import gzip
import os

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from camdweb.c2db.app import main
from camdweb.c2db.copy import copy_materials, main as copymain
from camdweb.test.c2db import create_tree
from camdweb.c2db.oqmd123 import db2json
from camdweb.test.html import check_html


@pytest.fixture
def oqmd_db_file(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp('tmp-oqmd')
    data = [('Mo', -11.202, 0.0),
            ('S48', -195.747, 0.0),
            ('Mo2S4', -44.212, -0.916)]
    path = tmp_path / 'oqmd123.db'
    db = connect(path)
    for formula, energy, hform in data:
        atoms = Atoms(formula)
        atoms.calc = SinglePointCalculator(atoms, energy=energy)
        db.write(atoms, hform=hform, uid=formula)
    gz = tmp_path / 'oqmd123.json.gz'
    db2json(path, gz)
    return gz


def test_c2db_missing_phonons(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    (tmp_path / 'MoS2/results-asr.phonons.json').unlink()
    (tmp_path / 'MoS2/results-asr.bader.json').unlink()

    copy_materials(tmp_path, ['MoS2*'], update_chull=False)
    # OK to do it again:
    copy_materials(tmp_path, ['MoS2*'], update_chull=False)

    with pytest.raises(FileNotFoundError):
        copy_materials(tmp_path, [])


def test_everything(oqmd_db_file):
    root = oqmd_db_file.parent
    create_tree(root)
    os.chdir(root)
    copymain([str(root), 'MoS2*'])

    app = main(['AB2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2/1'])
    assert len(app.materials) == 1

    html = app.index_page()
    check_html(html)

    # Compress one of the result files:
    bs = root / 'AB2/1MoS2/1/results-asr.bandstructure.json'
    with gzip.open(bs.with_suffix('.json.gz'), 'wt') as fd:
        fd.write(bs.read_text())
    bs.unlink()

    html = app.material_page('1MoS2-1')
    key = 'Charges [|e|]'
    passed = key in html
    assert passed
    check_html(html)

    (root / 'AB2/1MoS2/1/bader.json').unlink()
    (root / 'AB2/1MoS2/1/results-asr.shift.json').unlink()
    html = app.material_page('1MoS2-1')
    passed = key not in html  # Bader charge
    assert passed
