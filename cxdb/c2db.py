from __future__ import annotations

import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from ase import Atoms
from ase.io import read

from cxdb.asr_panel import ASRPanel, read_result_file
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.panel import Panel
from cxdb.web import CXDBApp
from cxdb.shift import ShiftPanel

RESULT_FILES = [
    'bandstructure',
    'phonons',
    'gs',
    'gs@calculate',
    'bader',
    'shift']


def copy_materials(path: Path, patterns: list[str]) -> None:
    names: defaultdict[str, int] = defaultdict(int)
    for pattern in patterns:
        print(pattern)
        for dir in path.glob(pattern):
            copy_material(dir, names)


def copy_material(dir: Path, names: defaultdict[str, int]) -> None:
    atoms = read(dir / 'gs.gpw')
    assert isinstance(atoms, Atoms)

    f = atoms.symbols.formula
    ab, xy, n = f.stoichiometry()
    name = f'{ab}/{n}{xy}'
    m = names[name] + 1
    names[name] = m
    folder = Path(name) / str(m)
    folder.mkdir(exist_ok=False, parents=True)

    atoms.write(folder / 'structure.xyz')

    # Copy result json-files:
    for name in RESULT_FILES:
        result = dir / f'results-asr.{name}.json'
        if result.is_file():
            shutil.copyfile(result, folder / result.name)

    def rrf(name: str) -> dict:
        return read_result_file(dir / f'results-asr.{name}.json')

    data = {}
    data['magstate'] = rrf('magstate')['magstate']
    data['has_inversion_symmetry'] = rrf(
        'structureinfo')['has_inversion_symmetry']
    gs = rrf('gs')
    data['gap'] = gs['gap']
    data['evac'] = gs['evac']
    data['hform'] = rrf('convex_hull')['hform']
    data['uid0'] = rrf('database.material_fingerprint')['uid']

    data['energy'] = atoms.get_potential_energy()

    (folder / 'data.json').write_text(json.dumps(data, indent=0))


class C2DBAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__(ndims=2)
        self.column_names.update(
            magstate='Magnetic state',
            ehull='Energy above convex hull [eV/atom]',
            hform='Heat of formation [eV/atom]',
            gap='Band gap (PBE) [eV]',
            energy='Energy [eV]',
            has_inversion_symmetry='Inversion symmetry',
            uid0='Old uid',
            evac='Vacuum level [eV]')
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            material.add_column(key, value)


def main(root: Path) -> CXDBApp:
    mlist: list[Material] = []
    for f in root.glob('A*/*/*/'):
        uid = f'{f.parent.name}-{f.name}'
        if len(mlist) % 20 == 0:
            print(end='.', flush=True)
        mlist.append(Material.from_file(f / 'structure.xyz', uid))
    print()

    panels: list[Panel] = [C2DBAtomsPanel()]
    for name in ['bandstructure',
                 'phonons',
                 'bader']:
        panels.append(ASRPanel(name))
    panels.append(ShiftPanel())

    materials = Materials(mlist, panels)

    initial_columns = ['magstate', 'ehull', 'hform', 'gap', 'formula', 'area']

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    if len(sys.argv) >= 3:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
    else:
        main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
