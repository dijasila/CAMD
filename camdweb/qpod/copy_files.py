"""Copy QPOD files from ASR-layout to uniform tree structure.

Tree structure::

   <host material>/<defects>/<charge>/
   <host material>/<pristine>/

Example::

   MoS2/v_Mo/charge_1/
   WS2/pristine_sc/
   ...

Build tree like this::

    $ cd /tmp
    $ mkdir tree
    $ cd tree
    $ python -m cxdb.qpod.copy_files <root-dir> <pattern> <pattern> ...

    """
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import rich.progress as progress
from camdweb.c2db.asr_panel import read_result_file
from camdweb.utils import read_atoms

RESULT_FILES = [
    'magstate',
    'magnetic_anisotropy',
    'structureinfo',
    # 'defect_symmetry', to be added; may or may not be present
    'gs',
    'gs@calculate',
    'database.material_fingerprint']

CHARGE_0_RESULT_FILES = ['defectinfo',
                         'sj_analyze']
# 'charge_neutrality'] handled in pris folder, same file
RESULT_FILES.extend(CHARGE_0_RESULT_FILES)


ROOT = Path('/home/niflheim2/cmr/defects/WIP/defects/rerun-qpod/tree-fabian2')

# Updated PATTERNS to match the new structure
PATTERNS = [
    'A*/*/*/defects.*/charge_*',
    'A*/*/*/defects.pristine*'
]


def copy_materials(root: Path, patterns: list[str]) -> None:
    dirs = [dir for pattern in patterns for dir in root.glob(pattern)]
    names: defaultdict[str, int] = defaultdict(int)

    with progress.Progress() as pb:
        pid = pb.add_task('Copying materials:', total=len(dirs))
        for dir in dirs:
            copy_material(dir, names)
            pb.advance(pid)


def copy_material(dir: Path, names: defaultdict[str, int]) -> None:
    gpw = dir / 'gs.gpw'
    if not gpw.is_file():
        return

    parts = str(dir).split('/')
    defects_index = next((i for i, part in enumerate(parts)
                          if 'defects.' in part), None)

    if defects_index is None or defects_index + 1 > len(parts):
        return

    # Handle pristine and defect folders
    material_folder = parts[defects_index - 1].split('-')[0]

    if 'defects.pristine' in str(dir):
        subfolder = parts[defects_index].split('.')[-1]
        folder = Path(material_folder) / subfolder
        RESULT_FILES.extend(['defectinfo'])  # 'charge_neutrality' conditional
        # charge_neutrality = dir / 'results-asr.charge_neutrality.json'
        RESULT_FILES.extend(['charge_neutrality'])
    else:
        defects_parts = parts[defects_index].split('.')
        subfolder = defects_parts[-1]
        charge_folder = parts[defects_index + 1]
        folder = Path(material_folder) / subfolder / charge_folder
        if charge_folder == 'charge_0':
            RESULT_FILES.extend(CHARGE_0_RESULT_FILES)

    def rrf(name: str) -> dict:
        return read_result_file(dir / f'results-asr.{name}.json')

    data = {}
    try:
        fingerprint = rrf('database.material_fingerprint')
        data['uid'] = fingerprint['uid']
    except FileNotFoundError:
        print("Database.fingerprint file not found", dir)
        return

    try:
        defectinfo = rrf('defectinfo')
        # Maybe take data from json directly at some later point
        data['host_name'] = defectinfo['host_name']
        data['defect_name'] = defectinfo['defect_name']
        data['charge_state'] = defectinfo['charge_state']
        data['host_gap_pbe'] = defectinfo['host_gap_pbe']
        data['host_gap_hse'] = defectinfo['host_gap_hse']
        data['host_hof'] = defectinfo['host_hof']
        data['host_uid'] = defectinfo['host_uid']
        data['host_spacegroup'] = defectinfo['host_spacegroup']
        data['host_pointgroup'] = defectinfo['host_pointgroup']
        # Defect-defect distance; 0 inplace of None:
        data['r_nn'] = defectinfo['R_nn'] or 0.0
    except FileNotFoundError:
        pass

    atoms = read_atoms(gpw)
    data['energy'] = atoms.get_potential_energy()

    folder.mkdir(exist_ok=False, parents=True)

    atoms.write(folder / 'structure.xyz')

    # Copy result json-files:
    for name in RESULT_FILES:
        result = dir / f'results-asr.{name}.json'

        if 'defects.pristine' in str(dir) and name == 'magstate':
            continue

        if result.is_file():
            shutil.copyfile(result, folder / result.name)

    (folder / 'data.json').write_text(json.dumps(data, indent=0))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        copy_materials(ROOT, PATTERNS)
    else:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
