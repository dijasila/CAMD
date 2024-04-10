import json
import sys
from collections import defaultdict
from pathlib import Path

from ase import Atoms
from ase.db import connect
from ase.formula import Formula


def copy_files(db_file: str) -> None:
    monolayers: defaultdict[str, dict[str, int]] = defaultdict(dict)
    for row in connect(db_file).select():
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
        if row.number_of_layers == 2:
            data['distance'] = distance(atoms)
            print(row.folder)
        (folder / 'data.json').write_text(json.dumps(data))
        for file in Path(row.folder).glob('results-asr.*.json'):
            (folder / file.name).write_bytes(file.read_bytes())


def distance(bilayer: Atoms) -> float:
    n = len(bilayer) // 2
    return bilayer.positions[n:, 2].min() - bilayer.positions[:n, 2].max()


if __name__ == '__main__':
    copy_files(sys.argv[1])
