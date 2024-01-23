import os

from cxdb.c2db.copy_files import copy_materials
from cxdb.c2db.app import main
from cxdb.test.c2db import create_tree


def test_c2db(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    copy_materials(tmp_path, ['MoS2*'])
    app = main(['AB2'])
    app.index()

    html = app.material('1MoS2-1')
    assert '1.24' in html  # Bader charge

    (tmp_path / 'AB2/1MoS2/1/results-asr.bader.json').unlink()
    (tmp_path / 'AB2/1MoS2/1/results-asr.shift.json').unlink()
    html = app.material('1MoS2-1')
    assert '1.24' not in html  # Bader charge


def test_c2db_missing_phonons(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    (tmp_path / 'MoS2/results-asr.phonons.json').unlink()
    copy_materials(tmp_path, ['MoS2*'])
