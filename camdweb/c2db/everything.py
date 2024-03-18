"""Create various tar-files and ASE db-files containing "everything".


"""
from ase.io import read
from pathlib import Path
import rich.progress as progress
from ase.db import connect


def create():
    db = connect('c2db.db')
    folders = list(Path().glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            uid = f'{f.parent.name}-{f.name}'
            atoms = read(f / 'structure.xyz')
            data = json.loads((f / 'data.json').read_text())

            pb.advance(pid)

    pool = None  # mp.Pool(maxtasksperchild=100)
