from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

from cxdb.asr_panel import ASRPanel
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp


def copy(pattern: str) -> None:
    names: defaultdict[str, int] = defaultdict(int)
    for dir in Path().glob(pattern):
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
                pass  # continue
            (folder / result.name).write_text(result.read_text())


def read_results(material, name):
    return json.loads(
        (material.folder / f'results-asr.{name}.json').read_text())


class C2DBAtomsPanel(AtomsPanel):
    column_names = AtomsPanel.column_names | {
        'area': 'Area [Å<sup>2</sup>]',
        'magstate': 'Magnetic',
        'ehull': 'Energy above convex hull [eV/atom]',
        'gap_pbe': 'Band gap (PBE) [eV]'}

    columns = list(column_names)

    def update_data(self, material):
        area = abs(np.linalg.det(material.atoms.cell[:2, :2]))
        material.add_column('area', area)
        magstate = read_results(material, 'magstate')['magstate']
        material.add_column('magstate', magstate)
        has_inversion_symmetry = read_results(
            material,
            'structureinfo')['kwargs']['data']['has_inversion_symmetry']
        material.add_column('has_inversion_symmetry', has_inversion_symmetry)


def main(root: Path) -> CXDBApp:
    mlist: list[Material] = []
    for f in root.glob('A*/*/*/'):
        uid = f'{f.parent.name}-{f.name}'
        if len(mlist) % 20 == 0:
            print(end='.', flush=True)
        mlist.append(Material(f, uid))
    print()

    panels = [C2DBAtomsPanel(),
              ASRPanel('bandstructure')]

    materials = Materials(mlist, panels)

    initial_columns = {'uid', 'energy', 'formula'}

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        copy(sys.argv[1])
    else:
        main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
