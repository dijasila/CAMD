"""Create various tar-files and ASE db-files containing "everything".

$ tar -czf c2db.tar.gz \
  --exclude "result*.json" --exclude "*.png" \
  C2DB/A* C2DB/convex-hulls/

"""
import json
from pathlib import Path

import rich.progress as progress
from ase.db import connect
from ase.io import read


def create():
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


if __name__ == '__main__':
    create()
