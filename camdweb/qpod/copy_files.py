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
    $ python -m camdweb.qpod.copy_files <root-dir> <pattern> <pattern> ...

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
    'gs',
    'gs@calculate']

CHARGE_0_RESULT_FILES = ['defectinfo',
                         'sj_analyze']
RESULT_FILES.extend(CHARGE_0_RESULT_FILES)


ROOT = Path('/home/niflheim2/cmr/WIP/defects/rerun-qpod')

PATTERNS = [
    'tree-*/A*/*/*/defects.*/charge_*',
    'tree-*/A*/*/*/defects.pristine*'
]

def copy_materials(root: Path, patterns: list[str]) -> None:
    dirs = [dir for pattern in patterns for dir in root.glob(pattern)]
    names: defaultdict[str, int] = defaultdict(int)
    
    assert root.exists()
    
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
        RESULT_FILES.extend(['defectinfo'])
        RESULT_FILES.extend(['charge_neutrality'])
    else:
        defects_parts = parts[defects_index].split('.')
        subfolder = defects_parts[-1]
        charge_folder = parts[defects_index + 1]
        folder = Path(material_folder) / subfolder / charge_folder
        RESULT_FILES.extend(['database.material_fingerprint'])
        
        if charge_folder == 'charge_0':
            RESULT_FILES.extend(CHARGE_0_RESULT_FILES)
        
        if (folder / 'results-asr.defect_symmetry.json').is_file():
            RESULT_FILES.extend(['defect_symmetry'])

    # def rrf(name: str) -> dict:
    #     return read_result_file(dir / f'results-asr.{name}.json')
    def rrf(name: str) -> dict:
        filepath = dir / f'results-asr.{name}.json'
        if not is_json_file_well_formatted(filepath):
            print(f"File {filepath} is not well formatted")
            if name == 'database.material_fingerprint':
                fix_json_file(filepath)
        
        return read_result_file(filepath)

    data = {}

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

    try:
        atoms = read_atoms(dir / 'structure.json')
    except FileNotFoundError:
        print('ERROR:', dir)

    if not 'defects.pristine' in str(dir):
        try:
            fingerprint = rrf('database.material_fingerprint')
            data['uid'] = fingerprint['uid']
        except FileNotFoundError:
            print("Database.fingerprint file not found", dir)
            return
        try:
            data['energy'] = atoms.get_potential_energy()
        except Exception:
            print(dir)

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

##### MAJOR HACKY FIX FOR *known* CURSED JSON FILES #####
##### This can be removed once the json files are fixed ##### 
def fix_json_file(filepath: str):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Check if the last line ends with '}}'
    if lines[-1].rstrip().endswith('}}'):
        # Remove the extra '}' from the last line
        lines[-1] = lines[-1].rstrip('}\n') + '\n'
        with open(filepath, 'w') as file:
            file.writelines(lines)
    
    try:
        with open(filepath, 'r') as file:
            json.load(file)
    except json.JSONDecodeError as e:
        error_message = str(e)
        if "Expecting ',' delimiter" in error_message:
            with open(filepath, 'a') as file:
                file.write('}')
        if "Extra data" in error_message:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            with open(filepath, 'w') as file:
                file.writelines(lines[:-1])

def is_json_file_well_formatted(filepath: str) -> bool:
    try:
        with open(filepath, 'r') as file:
            json.load(file)
        return True
    except json.JSONDecodeError:
        return False

##### END OF MAJOR HACKY FIX ##### 

if __name__ == '__main__':
    if len(sys.argv) == 1:
        copy_materials(ROOT, PATTERNS)
    else:
        copy_materials(Path(sys.argv[1]), sys.argv[2:])
