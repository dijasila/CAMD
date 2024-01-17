"""C2DB web-app.

This module has code to convert ~cmr/C2DB/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see cxdb.asr_panel module).
"""
from __future__ import annotations

import json
import shutil
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase import Atoms
from ase.io import read

from cxdb.asr_panel import ASRPanel, read_result_file
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.panel import Panel
from cxdb.shift import ShiftPanel
from cxdb.web import CXDBApp

RESULT_FILES = [
    'bandstructure',
    'phonons',
    'gs',
    'gs@calculate',
    'bader',
    'shift']

PATTERNS = [
    'tree/A*/*/*/',
    ]


def copy_all_c2db_materials():  # pragma: no cover
    """Copy C2DB files to uniform tree structure.

    Tree structure::

       <stoichiometry>/<formula-units><fomula>/<id>/

    Example::

       AB2/1MoS2/1/
       AB2/1MoS2/2/
       ...

    Build tree like this::

        $ cd /tmp
        $ mkdir tree
        $ cd tree
        $ python -c "from cxdb.c2db import *; copy_all_c2db_materials()"

    """
    #root = Path('/home/niflheim2/cmr/C2DB-ASR')
    root = Path('/home/tara/webpage/tree-asr')
    copy_materials(root, PATTERNS)


def copy_materials(root: Path, patterns: list[str]) -> None:
    dirs = [dir
            for pattern in patterns
            for dir in root.glob(pattern)
            if dir.name[0] != '.']
    names: defaultdict[str, int] = defaultdict(int)
    with progress.Progress() as pb:
        pid = pb.add_task('Copying matrerials:', total=len(dirs))
        for dir in dirs:
            copy_material(dir, names)
            pb.advance(pid)


def copy_material(dir: Path, names: defaultdict[str, int]) -> None:
    gpw = dir / 'gs.gpw'
    if not gpw.is_file():
        return  # pragma: no cover
    atoms = read(gpw)
    assert isinstance(atoms, Atoms)

    # Find uid (example: "1MoS2-1"):
    f = atoms.symbols.formula
    ab, xy, n = f.stoichiometry()
    name = f'{ab}/{n}{xy}'
    m = names[name] + 1
    names[name] = m
    folder = Path(name) / str(m)

    def rrf(name: str) -> dict:
        return read_result_file(dir / f'results-asr.{name}.json')

    data = {}
    try:
        data['magstate'] = rrf('magstate')['magstate']
        data['has_inversion_symmetry'] = rrf(
            'structureinfo')['has_inversion_symmetry']
        gs = rrf('gs')
        data['gap'] = gs['gap']
        data['evac'] = gs['evac']
        data['hform'] = rrf('convex_hull')['hform']
        data['uid0'] = rrf('database.material_fingerprint')['uid']
    except FileNotFoundError:  # pragma: no cover
        return

    data['energy'] = atoms.get_potential_energy()

    folder.mkdir(exist_ok=False, parents=True)

    atoms.write(folder / 'structure.xyz')

    # Copy result json-files:
    for name in RESULT_FILES:
        result = dir / f'results-asr.{name}.json'
        if result.is_file():
            shutil.copyfile(result, folder / result.name)

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
    """Create C2DB app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(files))
        for f in files:
            uid = f'{f.parent.name}-{f.name}'
            mlist.append(Material.from_file(f / 'structure.xyz', uid))
            pb.advance(pid)

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
    main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
