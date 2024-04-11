import os

from ase import Atoms
from ase.db import connect

from camdweb.bidb.app import main
from camdweb.bidb.copy import copy_files


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
                 binding_energy_gs=15.0,
                 folder=str(tmp_path))
        atoms = Atoms('H', pbc=(1, 1, 0))
        atoms.center(vacuum=1)
        db.write(atoms,
                 number_of_layers=1,
                 monolayer_uid='H-xyz',
                 folder=str(tmp_path))
    copy_files(dbfile)
    app = main(tmp_path)
    app.index_page()
    app.material_page('1H-1-stacking')
    html = app.material_page('1H-1')
    assert '15.000' in html
