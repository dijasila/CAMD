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
    $ python -m camdweb.c2db.copy <root-dir> <pattern> <pattern> ...

"""
from __future__ import annotations

import argparse
import json
import shutil
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from ase import Atoms
from ase.formula import Formula
from ase.io import read
from camdweb.c2db.emass import get_emass_data
import numpy as np

from camdweb import ColVal
from camdweb.c2db.asr_panel import read_result_file
from camdweb.c2db.convex_hull import update_chull_data
from camdweb.c2db.oqmd123 import read_oqmd123_data
from camdweb.utils import process_pool

RESULT_FILES = [
    'stiffness',
    'phonons',
    'deformationpotentials',
    'bandstructure',
    'pdos',
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
    'shift',
    'structureinfo',
    'collect_spiral',
    'dmi',
    'spinorbit',
    'spinorbit@calculate']

ROOT = Path('/home/niflheim2/cmr/C2DB-ASR')

PATTERNS = [
    'tree/A*/*/*/',
    'ICSD-COD/*el/*/',
    'adhoc_materials/*/',
    'tree_LDP/A*/*/*/',
    'tree_CDVAE/A*/*/*/',
    'tree_intercalated/A*/*/*/',
    'push-manti-tree/A*/*/*/',
    'tree_Wang23/A*/*/*/']


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.int64):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def copy_materials(root: Path,
                   patterns: list[str],
                   update_chull: bool = True,
                   processes: int = 1) -> None:
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
                print('No fingerprint:', fp)
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

    with process_pool(processes) as pool:
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
    """Used by Pool."""
    copy_material(*args)


def copy_material(fro: Path,
                  to: Path,
                  olduid: str,
                  uid: str) -> None:  # pragma: no cover
    structure_file = fro / 'structure.json'
    if not structure_file.is_file():
        return
    atoms = read(structure_file)
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
    data['has_inversion_symmetry'] = structure['has_inversion_symmetry']
    data['layergroup'] = structure.get('layergroup', '?')
    data['lgnum'] = structure.get('lgnum', -1)

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
        data['gap_dir'] = gs['gap_dir']
        data['gap_dir_nosoc'] = gs['gap_dir_nosoc']
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
            data[f'gap_{x}'] = r.get(f'gap_{x}', 0.0)
            data[f'gap_dir_{x}'] = r.get(f'gap_dir_{x}', 0.0)
            data[f'vbm_{x}'] = r.get(f'vbm_{x}')
            data[f'cbm_{x}'] = r.get(f'cbm_{x}')

    try:
        bse = rrf('bse')
    except FileNotFoundError:
        pass
    else:  # pragma: no cover
        data['E_B'] = bse.get('E_B')

    try:
        pol = rrf('polarizability')
    except FileNotFoundError:
        pass
    else:  # pragma: no cover
        for a in 'xyz':
            data[f'alpha{a}_el'] = pol[f'alpha{a}_el']

    try:
        irpol = rrf('infraredpolarizability')
    except FileNotFoundError:
        pass
    else:  # pragma: no cover
        for a in 'xyz':
            data[f'alpha{a}_lat'] = irpol.get(f'alpha{a}_lat')
            if f'alpha{a}_el' in data:
                data[f'alpha{a}'] = (
                    data[f'alpha{a}_el'] +  # type: ignore[operator]
                    data[f'alpha{a}_lat'])

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

    try:
        emass_data = rrf('effective_masses')
    except FileNotFoundError:
        pass
    else:
        emass_webpanel_data = get_emass_data(emass_data, atoms)
        with open(to / 'emass.json', 'w') as file:
            json.dump(emass_webpanel_data, file, indent=4, cls=NumpyEncoder)

    # Copy result json-files:
    for name in RESULT_FILES:
        result = fro / f'results-asr.{name}.json'
        target = to / result.name
        gzipped = target.with_suffix('.json.gz')
        if result.is_file() and not (target.is_file() or gzipped.is_file()):
            shutil.copyfile(result, target)

    path = to / 'data.json'
    if path.is_file():
        olddata = json.loads(path.read_text())
        data['ehull'] = olddata.get('ehull')
        data['hform'] = olddata.get('hform')

    # Remove None values:
    data = {key: value for key, value in data.items() if value is not None}

    path.write_text(json.dumps(data, indent=0))


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser()
    parser.add_argument('root', help='Root of ASR-tree to copy from.')
    patterns = ', '.join(f'"{p}"' for p in PATTERNS)
    parser.add_argument(
        'pattern', nargs='+',
        help='Glob pattern like "tree/A*/*/*/". '
        f'Use "ALL" to get all the standard patterns: {patterns}.')
    parser.add_argument('-s', '--skip-convex-hulls', action='store_true')
    parser.add_argument('-p', '--processes', type=int, default=1)
    args = parser.parse_args(argv)
    if args.pattern == ['ALL']:  # pragma: no cover
        args.pattern = PATTERNS
    copy_materials(Path(args.root), args.pattern, not args.skip_convex_hulls,
                   args.processes)


if __name__ == '__main__':
    main()
