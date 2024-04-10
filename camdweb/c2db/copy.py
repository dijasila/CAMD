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
from camdweb.c2db.emass import get_emass_data
import numpy as np

from camdweb import ColVal
from camdweb.c2db.asr_panel import read_result_file
from camdweb.c2db.convex_hull import update_chull_data
from camdweb.c2db.oqmd123 import read_oqmd123_data
from camdweb.utils import process_pool, read_atoms

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
    'raman',
    'bse',
    'piezoelectrictensor',
    'plasmafrequency',
    'gs',
    'gs@calculate',
    'fermisurface',
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


def all_dirs(root: Path,
             patterns: list[str]) -> list[Path]:
    return [dir
            for pattern in patterns
            for dir in root.glob(pattern)
            if dir.name[0] != '.']


def atoms_to_uid_name(atoms: Atoms) -> str:
    f = atoms.symbols.formula
    stoi, reduced, nunits = f.stoichiometry()
    return f'{nunits}{reduced}'


def create_uids(root: Path = ROOT,
                patterns: list[str] = PATTERNS) -> None:
    names: defaultdict[str, int] = defaultdict(int)
    new = []
    for dir in all_dirs(root, patterns):
        print(dir)
        try:
            atoms = read_atoms(dir / 'structure.json')
        except FileNotFoundError:
            print('ERROR:', dir)
            continue
        energy = atoms.get_potential_energy()
        uid = None
        olduid = None
        try:
            uid_data = json.loads((dir / 'uid.json').read_text())
        except FileNotFoundError:
            pass
        else:
            uid = uid_data['uid']
            olduid = uid_data.get('olduid')

        if olduid is None:
            fp = dir / 'results-asr.database.material_fingerprint.json'
            try:
                olduid = read_result_file(fp)['uid']
            except FileNotFoundError:
                pass

        name = atoms_to_uid_name(atoms)
        if uid:
            number = int(uid.split('-')[1])
            names[name] = max(names[name], number)
        else:
            new.append((name, energy, dir, olduid))

    print(len(new))
    for name, _, dir, olduid in sorted(new):
        number = names[name] + 1
        uid = f'{name}-{number}'
        names[name] = number
        uid_data = {'uid': uid}
        if olduid:
            uid_data['olduid'] = olduid
        (dir / 'uid.json').write_text(json.dumps(uid_data, indent=2))


def copy_materials(root: Path,
                   patterns: list[str],
                   update_chull: bool = True,
                   processes: int = 1) -> None:
    dirs = all_dirs(root, patterns)
    print(len(dirs), 'folders')

    names: defaultdict[str, int] = defaultdict(int)
    work = []
    with progress.Progress() as pb:
        pid = pb.add_task('Finding UIDs:', total=len(dirs))
        for dir in dirs:
            try:
                uid = json.loads((dir / 'uid.json').read_text())['uid']
            except FileNotFoundError:
                try:
                    name = atoms_to_uid_name(
                        read_atoms(dir / 'structure.json'))
                except FileNotFoundError:
                    pb.advance(pid)
                    continue
                names[name] += 1
                number = names[name]
                uid = f'{name}-{number}t'
            name, tag = uid.split('-')
            stoichiometry, _, nunits = Formula(name).stoichiometry()
            folder = Path(f'{stoichiometry}/{name}/{tag}')
            work.append((dir, folder, uid))
            pb.advance(pid)

    parent_folders = set()
    for _, folder, _ in work:
        parent_folders.add(folder.parent)
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

    # Put logo in the right place:
    logo = Path('c2db-logo.png')
    if not logo.is_file():
        shutil.copyfile(Path(__file__).parent / 'logo.png', logo)


def worker(args):  # pragma: no cover
    """Used by Pool."""
    copy_material(*args)


def copy_material(fro: Path,
                  to: Path,
                  uid: str) -> None:  # pragma: no cover
    structure_file = fro / 'structure.json'
    if not structure_file.is_file():
        return
    atoms = read_atoms(structure_file)

    to.mkdir(exist_ok=True)

    def rrf(name: str) -> dict:
        return read_result_file(fro / f'results-asr.{name}.json')

    # None values will be removed later
    data: dict[str, ColVal | None] = {'folder': str(fro)}

    try:
        data['magstate'] = rrf('magstate')['magstate']
    except FileNotFoundError:
        pass

    try:
        data['spin_axis'] = rrf('magnetic_anisotropy')['spin_axis']
    except FileNotFoundError:
        pass

    # Read uid and perhaps olduid:
    try:
        data.update(json.loads((fro / 'uid.json').read_text()))
    except FileNotFoundError:
        data['uid'] = uid
    else:
        assert data['uid'] == uid

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
        maxphononfreq = ph['omega_kl'][0].max()

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
    else:
        for a in 'xyz':
            data[f'alpha{a}_el'] = pol[f'alpha{a}_el']
        dct = {'omega_w': pol['frequencies'].tolist()}
        for v in 'xyz':
            alpha = pol[f'alpha{v}_w']
            dct[f'alpha{v}_re_w'] = alpha.real.tolist()
            dct[f'alpha{v}_im_w'] = alpha.imag.tolist()
        (to / 'opt-polarizability.json').write_text(json.dumps(dct))

    try:
        irpol = rrf('infraredpolarizability')
    except FileNotFoundError:
        pass
    else:
        for a in 'xyz':
            data[f'alpha{a}_lat'] = irpol.get(f'alpha{a}_lat')
            if f'alpha{a}_el' in data:
                data[f'alpha{a}'] = (
                    data[f'alpha{a}_el'] +  # type: ignore[operator]
                    data[f'alpha{a}_lat'])
        alpha_wvv = irpol['alpha_wvv']
        dct = {'maxphononfreq': maxphononfreq,
               'omega_w': irpol['omega_w'].tolist(),
               'alpha_re_wvv': alpha_wvv.real.tolist(),
               'alpha_im_wvv': alpha_wvv.imag.tolist()}
        (to / 'ir-polarizability.json').write_text(json.dumps(dct))

    try:
        dct = rrf('plasmafrequency')
    except FileNotFoundError:
        pass
    else:
        for v in 'xy':
            data[f'plasmafrequency_{v}'] = dct[f'plasmafrequency_{v}']

    data['energy'] = atoms.get_potential_energy()

    try:
        info = json.loads((fro / 'info.json').read_text())
    except FileNotFoundError:
        pass
    else:
        for key in ['icsd_id', 'cod_id', 'doi']:
            if key in info:
                data[key] = info[key]

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
    parser.add_argument('root', help='Root of ASR-tree to copy from. '
                        'Example: "~cmr/C2DB-ASR/".')
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
