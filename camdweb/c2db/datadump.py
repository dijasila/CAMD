"""Create various tar-files and ASE db-files containing "everything".

$ tar -czf c2db.tar.gz \
  --exclude "result*.json" --exclude "*.png" \
  C2DB/A* C2DB/convex-hulls/

"""
import json
from pathlib import Path
import shutil

import rich.progress as progress
from ase.db import connect
from ase.io import read


def create_db_file():
    db = connect('c2db.db')
    folders = list(Path().glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            # uid = f'{f.parent.name}-{f.name}'
            atoms = read(f / 'structure.xyz')
            data = json.loads((f / 'data.json').read_text())
            data.pop('energy')
            db.write(atoms, **data)
            pb.advance(pid)


def create_data_dump():
    folders = list(Path().glob('A*/*/*/'))
    data = Path('data')
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            to = data / f
            to.mkdir(parents=True, exist_ok=True)
            for name in ['structure.xyz',
                         'data.json',
                         'results-asr.phonons.json']:
                fro = f / name
                if fro.is_file():
                    shutil.copy(fro, to)
            pb.advance(pid)


if __name__ == '__main__':
    # create_db_file()
    create_data_dump()

