from __future__ import annotations
import argparse
from pathlib import Path

from ase import Atoms
from ase.io import read

from camdweb.material import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.web import CAMDApp


class MyAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()
        self.column_names.update(
            energy='Energy [eV]',
            fmax='Maximum force [eV/Å]',
            smax='Maximum stress component [eV/Å<sup>3</sup>]',
            magmom='Total magnetic moment [μ<sub>B</sub>]')
        self.columns += ['energy', 'fmax', 'smax', 'magmom']

    def update_data(self, material):
        super().update_data(material)
        atoms = material.atoms
        try:
            results = atoms.calc.results
        except AttributeError:
            return
        energy = results.get('energy')
        if energy is not None:
            material.add_column('energy', energy)
        forces = results.get('forces')
        if forces is not None:
            material.add_column('fmax', (forces**2).sum(axis=1).max()**0.5)
        stress = results.get('stress')
        if stress is not None:
            material.add_column('smax', abs(stress).max())
        magmom = results.get('magmom')
        if magmom is not None:
            material.add_column('magmom', magmom)


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
        rows.append(Material(path.parent, str(i), atoms))

    panels = [MyAtomsPanel()]
    materials = Materials(rows, panels)

    initial_columns = ['uid', 'formula', 'energy', 'fmax', 'smax', 'magmom']
    for key in ['length', 'area', 'volume']:
        if key in materials.column_names:
            initial_columns.append(key)

    root = Path.cwd()

    app = CAMDApp(materials, initial_columns, root)
    if run:  # pragma: no cover
        app.app.run(host='0.0.0.0', port=8080, debug=True)
    return app


if __name__ == '__main__':
    main()
