import os
from pathlib import Path

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect

from cxdb.c2db.app import main
from cxdb.c2db.convex_hull import update_chull_data
from cxdb.c2db.copy_files import copy_materials
from cxdb.test.c2db import create_tree
from cxdb.c2db.oqmd123 import main as create_oqmd_json_gz_file


def create_oqmd_db_file(dir: Path):
    data = [('Mo', -11.202, 0.0),
            ('S48', -195.747, 0.0),
            ('Mo2S4', -44.212, -0.916)]
    db = connect(dir / 'oqmd123.db')
    for formula, energy, hform in data:
        atoms = Atoms(formula)
        atoms.calc = SinglePointCalculator(atoms, energy=energy)
        db.write(atoms, hform=hform, uid=formula)


def test_c2db_missing_phonons(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    (tmp_path / 'MoS2/results-asr.phonons.json').unlink()
    copy_materials(tmp_path, ['MoS2*'])


def test_everything(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    copy_materials(tmp_path, ['MoS2*'])
    create_oqmd_db_file(tmp_path)
    create_oqmd_json_gz_file(tmp_path / 'oqmd123.db')
    update_chull_data(tmp_path)

    app = main(['AB2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2/1'])
    assert len(app.materials) == 1

    app.index()

    with gzip.open('...'):

    html = app.material('1MoS2-1')
    assert '1.24' in html  # Bader charge

    (tmp_path / 'AB2/1MoS2/1/results-asr.bader.json').unlink()
    (tmp_path / 'AB2/1MoS2/1/results-asr.shift.json').unlink()
    (tmp_path / 'AB2/1MoS2/1/results-asr.convex_hull.json').unlink()
    html = app.material('1MoS2-1')
    assert '1.24' not in html  # Bader charge
