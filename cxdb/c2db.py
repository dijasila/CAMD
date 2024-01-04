from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

from ase import Atoms
from ase.io import read

from cxdb.asr_panel import ASRPanel
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp


def copy(path: Path, pattern: str) -> None:
    names: defaultdict[str, int] = defaultdict(int)
    print(pattern)
    for dir in path.glob(pattern):
        atoms = read(dir / 'gs.gpw')
        assert isinstance(atoms, Atoms)
        f = atoms.symbols.formula
        ab, xy, n = f.stoichiometry()
        name = f'{ab}/{n}{xy}'
        m = names[name] + 1
        names[name] = m
        folder = Path(name) / str(m)
        folder.mkdir(exist_ok=True, parents=True)
        atoms.write(folder / 'structure.xyz')
        for result in dir.glob('results-asr.*.json'):
            if '@' in result.name:
                if result.name != 'results-asr.gs@calculate.json':
                    continue
            (folder / result.name).write_text(result.read_text())


def read_results(material, name, hmm=False):
    dct = json.loads(
        (material.folder / f'results-asr.{name}.json').read_text())
    if hmm:
        dct = dct['kwargs']['data']
    return dct


class C2DBAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__(2)
        self.column_names.update(
            magstate='Magnetic',
            ehull='Energy above convex hull [eV/atom]',
            hform='Heat of formation [eV/atom]',
            gap='Band gap (PBE) [eV]',
            energy='Energy [eV]',
            has_inversion_symmetry='Inversion symmetry')
        self.columns = list(self.column_names)

    def update_data(self, material):
        super().update_data(material)
        energy = material.atoms.get_potential_energy()
        material.add_column('energy', energy)
        magstate = read_results(material, 'magstate')['magstate']
        material.add_column('magstate', magstate)
        has_inversion_symmetry = read_results(
            material, 'structureinfo', 1)['has_inversion_symmetry']
        material.add_column('has_inversion_symmetry', has_inversion_symmetry)
        gap = read_results(material, 'gs', 1)['gap']
        material.add_column('gap', gap)
        hform = read_results(material, 'convex_hull', 1)['hform']
        material.add_column('hform', hform)


def main(root: Path) -> CXDBApp:
    mlist: list[Material] = []
    for f in root.glob('A*/*/*/'):
        uid = f'{f.parent.name}-{f.name}'
        if len(mlist) % 20 == 0:
            print(end='.', flush=True)
        mlist.append(Material(f, uid))
    print()

    panels = [C2DBAtomsPanel(),
              ASRPanel('bandstructure'),
              ASRPanel('phonons')]

    materials = Materials(mlist, panels)

    initial_columns = {'magstate', 'ehull', 'hform', 'gap', 'formula', 'area'}

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        copy(Path(sys.argv[1]), sys.argv[2])
    else:
        main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
