from cxdb.bidb import main, expand
from ase.db import connect
from ase import Atoms
import os


def test_bidb(tmp_path):
    os.chdir(tmp_path)
    dbfile = tmp_path / 'bidb.db'
    with connect(dbfile) as db:
        atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=(1, 1, 0))
        atoms.center(vacuum=1)
        db.write(atoms,
                 number_of_layers=2,
                 monolayer_uid='H-xyz',
                 bilayer_uid='H-xyz-stacking',
                 cod_id='A23462346',
                 extra=27,
                 binding_energy_zscan=15.0)
        atoms = Atoms('H', pbc=(1, 1, 0))
        atoms.center(vacuum=1)
        db.write(atoms,
                 number_of_layers=1,
                 monolayer_uid='H-xyz')
    expand(dbfile)
    app = main(tmp_path)
    app.index()
    print(app.materials._materials)
    app.material('1H-0-stacking')
    html = app.material('1H-0')
    assert '15.000' in html
