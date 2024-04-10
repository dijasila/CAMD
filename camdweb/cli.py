from __future__ import annotations
import argparse
from pathlib import Path

from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.web import CAMDApp
from camdweb.html import table
from camdweb.utils import read_atoms

COLUMN_DESCRIPTIONS = dict(
    energy='Energy [eV]',
    fmax='Maximum force [eV/Å]',
    smax='Maximum stress component [eV/Å<sup>3</sup>]',
    magmom='Total magnetic moment [μ<sub>B</sub>]')


class MyAtomsPanel(AtomsPanel):
    column_descriptions = COLUMN_DESCRIPTIONS

    def update_material(self, material: Material) -> None:
        try:
            results = material.atoms.calc.results
        except AttributeError:
            return
        energy = results.get('energy')
        if energy is not None:
            material.columns['energy'] = energy
        forces = results.get('forces')
        if forces is not None:
            material.columns['fmax'] = (forces**2).sum(axis=1).max()**0.5
        stress = results.get('stress')
        if stress is not None:
            material.columns['smax'] = abs(stress).max()
        magmom = results.get('magmom')
        if magmom is not None:
            material.columns['magmom'] = magmom

    def create_column_one(self,
                          material: Material) -> str:
        rows = self.table_rows(
            material,
            ['formula', 'energy', 'fmax', 'smax', 'magmom',
             'length', 'area', 'volume'])
        return table(None, rows)


def main(argv: list[str] | None = None,
         run=True) -> CAMDApp:
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='+',
                        help='Filename of atomic structure file.')
    args = parser.parse_args(argv)
    rows: list[Material] = []
    for i, filename in enumerate(args.filename):
        path = Path(filename)
        atoms = read_atoms(path)
        rows.append(Material(str(i), atoms))

    panels = [MyAtomsPanel()]
    materials = Materials(rows, panels)

    initial_columns = ['uid', 'formula', 'energy', 'fmax', 'smax', 'magmom']
    for key in ['length', 'area', 'volume']:
        if key in materials.column_descriptions:
            initial_columns.append(key)

    root = Path.cwd()

    app = CAMDApp(materials, initial_columns, root=root)
    if run:  # pragma: no cover
        app.app.run(host='0.0.0.0', port=8080, debug=True)
    return app


if __name__ == '__main__':
    main()
