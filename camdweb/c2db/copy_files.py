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
    $ python -m camdweb.c2db.copy_files <root-dir> <pattern> <pattern> ...

    """
import json
import multiprocessing as mp
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase import Atoms
from ase.formula import Formula
from ase.io import read

from camdweb import ColVal
from camdweb.c2db.asr_panel import read_result_file
from camdweb.c2db.convex_hull import update_chull_data
from camdweb.c2db.oqmd123 import read_oqmd123_data

RESULT_FILES = [
    'stiffness',
    'phonons',
    'deformationpotentials',
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


def copy_materials(root: Path, patterns: list[str],
                   update_chull: bool = True) -> None:
    try:
        uids = json.loads(Path('uids.json').read_text())
    except FileNotFoundError:
        uids = {}
    print(len(uids), 'UIDs')

    names: defaultdict[str, int] = defaultdict(int)
    for uid in uids.values():
        name = uid.split('-')[0]
        names[name] += 1

    dirs = [dir
            for pattern in patterns
            for dir in root.glob(pattern)
            if dir.name[0] != '.']
    print(len(dirs), 'folders')

    work = []
    parent_folders = set()
    with progress.Progress() as pb:
        pid = pb.add_task('Finding UIDs:', total=len(dirs))
        for dir in dirs:
            fp = dir / 'results-asr.database.material_fingerprint.json'
            try:
                olduid = read_result_file(fp)['uid']
            except FileNotFoundError:  # pragma: no cover
                print(fp)
                continue
            f = Formula(olduid.split('-')[0])
            stoi, reduced, nunits = f.stoichiometry()
            name = f'{nunits}{reduced}'
            uid = uids.get(olduid)
            if uid is not None:
                name1, x = uid.split('-')
                number = int(x)
                assert name1 == name
            else:
                names[name] += 1
                number = names[name]
                uid = f'{nunits}{reduced}-{number}'
                uids[olduid] = uid
            folder = Path(f'{stoi}/{name}/{number}')
            parent_folders.add(folder.parent)
            work.append((dir, folder, olduid, uid))
            pb.advance(pid)

    Path('uids.json').write_text(json.dumps(uids, indent=1))

    for folder in parent_folders:
        folder.mkdir(exist_ok=True, parents=True)

    with mp.Pool(processes=8) as pool:
        print(pool)
        with progress.Progress() as pb:
            pid = pb.add_task('Copying materials:', total=len(work))
            for _ in pool.imap_unordered(worker, work):
                pb.advance(pid)

    if update_chull:
        # Calculate hform, ehull, ...
        try:
            oqmd_path = Path('oqmd123.json.gz')
            atomic_energies, refs = read_oqmd123_data(oqmd_path)
        except FileNotFoundError:
            raise FileNotFoundError(
                f'Could not find {oqmd_path}.\n'
                'Please download the oqmd123.db file:\n\n'
                '   wget https://cmr.fysik.dtu.dk/_downloads/oqmd123.db\n\n'
                'and convert to json with:\n\n'
                '   python -m camdweb.c2db.oqmd123 <path-to-oqmd123.db>\n')
        update_chull_data(atomic_energies, refs)


def worker(args):  # pragma: no cover
    """Used by Pool"""
    copy_material(*args)


def copy_material(fro: Path,
                  to: Path,
                  olduid: str,
                  uid: str) -> None:  # pragma: no cover
    gpw = fro / 'gs.gpw'
    if gpw.is_file():
        atoms = read(gpw)
    else:
        atoms = read(fro / 'structure.json')
    assert isinstance(atoms, Atoms)

    def rrf(name: str) -> dict:
        return read_result_file(fro / f'results-asr.{name}.json')

    # None values will be removed later:
    data: dict[str, ColVal | None] = {
        'uid': uid,
        'olduid': olduid,
        'folder': str(fro)}
    try:
        data['magstate'] = rrf('magstate')['magstate']
    except FileNotFoundError:
        pass
    try:
        data['spin_axis'] = rrf('magnetic_anisotropy')['spin_axis']
    except FileNotFoundError:
        pass

    structure = rrf('structureinfo')
    for key in ['has_inversion_symmetry', 'layergroup', 'lgnum']:
        data[key] = structure[key]

    try:
        data['label'] = rrf('c2db.labels')['label']
    except FileNotFoundError:
        pass

    try:
        gs = rrf('gs')
    except FileNotFoundError:
        pass
    else:
        data['gap'] = gs['gap']
        data['evac'] = gs['evac']
        data['efermi'] = gs['efermi']

    try:
        ph = rrf('phonons')
    except FileNotFoundError:
        dyn_stab_phonons = 'unknown'
    else:
        data['minhessianeig'] = ph['minhessianeig']
        dyn_stab_phonons = ph['dynamic_stability_phonons']

    try:
        ph = rrf('stiffness')
    except FileNotFoundError:
        dyn_stab_stiffness = 'unknown'
    else:  # pragma: no cover
        dyn_stab_stiffness = ph['dynamic_stability_stiffness']

    data['dyn_stab'] = (dyn_stab_phonons == 'high' and
                        dyn_stab_stiffness == 'high')

    for x in ['hse', 'gw']:
        try:
            r = rrf(x)
        except FileNotFoundError:
            pass
        else:  # pragma: no cover
            data[f'gap_{x}'] = r.get(f'gap_{x}')
            data[f'gap_dir_{x}'] = r.get(f'gap_dir_{x}')
            data[f'vbm_{x}'] = r.get(f'vbm_{x}')
            data[f'cbm_{x}'] = r.get(f'cbm_{x}')

    try:
        pol = rrf('polarizability')
    except FileNotFoundError:
        pass
    else:  # pragma: no cover
        for a in 'xyz':
            data[f'alpha{a}_el'] = pol[f'alpha{a}_el']
            data[f'alpha{a}_lat'] = pol.get(f'alpha{a}_lat')

    data['energy'] = atoms.get_potential_energy()

    try:
        info = json.loads((fro / 'info.json').read_text())
    except FileNotFoundError:
        pass
    else:
        for key in ['icsd_id', 'cod_id', 'doi']:
            if key in info:
                data[key] = info[key]

    to.mkdir(exist_ok=True)

    atoms.write(to / 'structure.xyz')

    try:
        bc = rrf('bader')['bader_charges']
    except FileNotFoundError:
        pass
    else:
        (to / 'bader.json').write_text(
            json.dumps({'charges': bc.tolist()}))

    # Copy result json-files:
    for name in RESULT_FILES:
        result = fro / f'results-asr.{name}.json'
        if result.is_file():
            shutil.copyfile(result, to / result.name)

    data = {key: value for key, value in data.items() if value is not None}
    (to / 'data.json').write_text(json.dumps(data, indent=0))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        copy_materials(ROOT, PATTERNS)
    else:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
