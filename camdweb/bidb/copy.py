import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from ase import Atoms
from ase.db import connect
from ase.formula import Formula


def copy_files(db_file: str,
               folders_file: str,
               home: str = '/home/') -> None:
    folder_dict = json.loads(Path(folders_file).read_text())
    include_folders = set()
    for monolayer_folder, bilayer_folders in folder_dict.items():
        include_folders.add(monolayer_folder)
        include_folders.update(bilayer_folders)
    monolayers: defaultdict[str, dict[str, int]] = defaultdict(dict)
    for row in connect(db_file).select():
        if row.folder not in include_folders:
            continue
        f = Formula(row.formula)
        ab, xy, n = f.stoichiometry()
        n //= row.number_of_layers
        name = f'{ab}/{n}{xy}'
        ids = monolayers[name]
        if row.monolayer_uid in ids:
            i = ids[row.monolayer_uid]
        else:
            i = len(ids) + 1
            ids[row.monolayer_uid] = i
        folder = Path(name + f'-{i}')
        if row.number_of_layers == 1:
            folder /= 'monolayer'
        else:
            folder /= row.bilayer_uid.split('-', 2)[2]
        folder.mkdir(exist_ok=True, parents=True)
        atoms = row.toatoms()
        atoms.write(folder / 'structure.xyz')
        data = dict(row.key_value_pairs)
        for key in ['binding_energy_gs', 'binding_energy_zscan']:
            if key in data:
                if data[key] == '-':
                    del data[key]
                else:
                    data[key] *= 1000
        if row.number_of_layers == 2:
            data['distance'] = distance(atoms)
        (folder / 'data.json').write_text(json.dumps(data, indent=2))
        dir = Path(row.folder.replace('/home/', home))
        if 0:
            for file in dir.glob('results-asr.*.json'):
                (folder / file.name).write_bytes(file.read_bytes())
        for name in ['pdos']:
            file = dir / f'results-asr.{name}.json'
            if file.is_file():
                txt = file.read_text()
                try:
                    json.loads(txt)
                except json.JSONDecodeError:
                    print(file)
                else:
                    (folder / file.name).write_text(txt)

    # Put logo in the right place:
    logo = Path('bidb-logo.png')
    if not logo.is_file():
        shutil.copyfile(Path(__file__).parent / 'logo.png', logo)


def distance(bilayer: Atoms) -> float:
    n = len(bilayer) // 2
    return bilayer.positions[n:, 2].min() - bilayer.positions[:n, 2].max()


if __name__ == '__main__':
    copy_files(*sys.argv[1:])
