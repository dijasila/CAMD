from cxdb.cmr import main
from ase.db import connect
from ase import Atoms
import os


def test_cmr(tmp_path):
    os.chdir(tmp_path)
    dbfile = tmp_path / 'abc.db'
    with connect(dbfile) as db:
        atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
        atoms.center(vacuum=1)
        db.write(atoms, abc=27.3)
        atoms = Atoms('H')
        atoms.center(vacuum=1)
        db.write(atoms, abc=1.2)
    app = main([dbfile])
    app.overview()
    app.index('abc')
    app.material('abc', '1')
