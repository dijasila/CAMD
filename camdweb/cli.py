from __future__ import annotations
import argparse
from pathlib import Path

from ase import Atoms
from ase.io import read

from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.web import CAMDApp
from camdweb.html import table

COLUMN_DESCRIPTIONS = dict(
    energy='Energy [eV]',
    fmax='Maximum force [eV/Å]',
    smax='Maximum stress component [eV/Å<sup>3</sup>]',
    magmom='Total magnetic moment [μ<sub>B</sub>]')


class MyMaterial(Material):
    def __init__(self, folder, uid, atoms):
        super().__init__(folder, uid, atoms)
        try:
            results = atoms.calc.results
        except AttributeError:
            results = {}
        self.energy = results.get('energy')
        forces = results.get('forces')
        self.fmax = (None if forces is None else
                     (forces**2).sum(axis=1).max()**0.5)
        stress = results.get('stress')
        self.smax = None if stress is None else abs(stress).max()
        self.magmom = results.get('magmom')

    def get_columns(self):
        columns = super().get_columns()
        for key in ['energy', 'fmax', 'smax', 'magmom']:
            value = getattr(self, key, None)
            if value is not None:
                columns[key] = value
        return columns


class MyAtomsPanel(AtomsPanel):
    def create_column_one(self,
                          material: Material) -> str:
        rows = []
        for key in ['formula', 'energy', 'fmax', 'smax', 'magmom',
                    'length', 'area', 'volume']:
            value = getattr(material, key, None)
            if value is not None:
                desc = COLUMN_DESCRIPTIONS.get(key)
                desc = desc or COMMON_COLUMN_DESCRIPTIONS.get(key, key)
                rows.append([desc, material.html_format_column(key, value)])
        print(rows)
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
        atoms = read(path)
        assert isinstance(atoms, Atoms)
        rows.append(MyMaterial(path.parent, str(i), atoms))

    panels = [MyAtomsPanel()]
    materials = Materials(rows, panels, COLUMN_DESCRIPTIONS)

    initial_columns = ['uid', 'formula', 'energy', 'fmax', 'smax', 'magmom']
    for key in ['length', 'area', 'volume']:
        if key in materials.column_descriptions:
            initial_columns.append(key)

    root = Path.cwd()

    app = CAMDApp(materials, initial_columns, root)
    if run:  # pragma: no cover
        app.app.run(host='0.0.0.0', port=8080, debug=True)
    return app


if __name__ == '__main__':
    main()
