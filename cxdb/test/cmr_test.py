from cxdb.cmr.app import main
from ase.db import connect
from ase import Atoms
import os


def test_cmr(tmp_path):
    os.chdir(tmp_path)
    dbfile1 = tmp_path / 'abc.db'
    with connect(dbfile1) as db:
        atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
        atoms.center(vacuum=1)
        db.write(atoms, abc=27.3)
        atoms = Atoms('H')
        atoms.center(vacuum=1)
        db.write(atoms, abc=1.2)
    dbfile2 = tmp_path / 'solar.db'
    with connect(dbfile2) as db:
        atoms = Atoms('N2', [(0, 0, 0), (1.1, 0, 0)])
        atoms.center(vacuum=1)
        db.write(atoms, KS_gap=27.3, DeltaU=1.3)
    app = main([dbfile1, dbfile2])
    app.overview()
    app.index('abc')
    app.material('abc', '1')
    app.index('solar')
    app.material('solar', '1')
    app.download_db_file('abc')
    app.favicon()
    dct = app.callback('abc', dict(name='atoms', uid='1', data='1'))
    assert 'data' in dct
