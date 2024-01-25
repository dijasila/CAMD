"""Copy C2DB files from ASR-layout to uniform tree structure.

Tree structure::

   <stoichiometry>/<formula-units><formula>/<id>/

Example::

   AB2/1MoS2/1/
   AB2/1MoS2/2/
   ...

Build tree like this::

    $ cd /tmp
    $ mkdir tree
    $ cd tree
    $ python -m cxdb.c2db.copy_files <root-dir> <pattern> <pattern> ...

    """
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase import Atoms
from ase.io import read
from cxdb.c2db.asr_panel import read_result_file

RESULT_FILES = [
    'convex_hull',
    'stiffness',
    'phonons',
    # 'deformationpotentials',
    'bandstructure',
    'pdos',
    'effective_masses',
    'hse',
    'gw',
    'borncharges',
    'shg',
    'polarizability',
    'infraredpolarizability',
    'raman',
    'bse',
    'bader',
    'piezoelectrictensor',
    'gs',
    'gs@calculate',
    'shift']

ROOT = Path('/home/niflheim2/cmr/C2DB-ASR')

PATTERNS = [
    'tree/A*/*/*/',
    'ICSD-COD/*el/*/',
    'adhoc_materials/*/',
    'tree_LDP/A*/*/*/',
    'tree_CDVAE/A*/*/*/',
    'tree_intercalated/A*/*/*/',
    'push-manti-tree/A*/*/*/']
# '/home/niflheim2/pmely/trees_to_collect/tree_Wang23/A*/*/*/'


def copy_materials(root: Path, patterns: list[str]) -> None:
    dirs = [dir
            for pattern in patterns
            for dir in root.glob(pattern)
            if dir.name[0] != '.']
    names: defaultdict[str, int] = defaultdict(int)
    print(len(dirs), 'folders')
    with progress.Progress() as pb:
        pid = pb.add_task('Copying materials:', total=len(dirs))
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
    uid = f'{n}{xy}-{m}'

    def rrf(name: str) -> dict:
        return read_result_file(dir / f'results-asr.{name}.json')

    data = {'uid': uid}
    try:
        data['magstate'] = rrf('magstate')['magstate']
        data['has_inversion_symmetry'] = rrf(
            'structureinfo')['has_inversion_symmetry']
        gs = rrf('gs')
        data['gap'] = gs['gap']
        data['evac'] = gs['evac']
        data['uid0'] = rrf('database.material_fingerprint')['uid']
    except FileNotFoundError:  # pragma: no cover
        return

    try:
        data['minhessianeig'] = rrf('phonons')['minhessianeig']
    except FileNotFoundError:
        pass

    data['energy'] = atoms.get_potential_energy()

    folder.mkdir(exist_ok=False, parents=True)

    atoms.write(folder / 'structure.xyz')

    # Copy result json-files:
    for name in RESULT_FILES:
        result = dir / f'results-asr.{name}.json'
        if result.is_file():
            shutil.copyfile(result, folder / result.name)
    (folder / 'data.json').write_text(json.dumps(data, indent=0))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        copy_materials(ROOT, PATTERNS)
    else:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
