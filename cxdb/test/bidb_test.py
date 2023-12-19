from cxdb.bidb import main, expand
from ase.db import connect
from ase import Atoms
import os


def test_bidb(tmp_path):
    os.chdir(tmp_path)
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
    atoms.center(vacuum=1)
    dbfile = tmp_path / 'bidb.db'
    with connect(dbfile) as db:
        db.write(atoms,
                 number_of_layers=2,
                 monolayer_uid='H-xyz',
                 bilayer_uid='H-xyz-stacking')
    expand(dbfile)
    app = main(tmp_path)
    _ = app.index()
    _ = app.material('1H-0-stacking')
