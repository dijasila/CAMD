import gzip
import os

import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from camdweb.c2db.app import main
from camdweb.c2db.convex_hull import update_chull_data, read_chull_data
from camdweb.c2db.copy_files import copy_materials
from camdweb.test.c2db import create_tree


@pytest.fixture(scope='module')
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
    return path


def test_c2db_missing_phonons(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    (tmp_path / 'MoS2/results-asr.phonons.json').unlink()
    (tmp_path / 'MoS2/results-asr.bader.json').unlink()
    copy_materials(tmp_path, ['MoS2*'])


def test_everything(tmp_path, oqmd_db_file):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    copy_materials(tmp_path, ['MoS2*'])
    atomic_energies, refs = read_chull_data(oqmd_db_file.parent)
    update_chull_data(atomic_energies, refs, tmp_path)
    app = main(['AB2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2'])
    assert len(app.materials) == 1
    app = main(['AB2/1MoS2/1'])
    assert len(app.materials) == 1

    app.index()

    # Compress one of the result files:
    bs = tmp_path / 'AB2/1MoS2/1/results-asr.bandstructure.json'
    with gzip.open(bs.with_suffix('.json.gz'), 'wt') as fd:
        fd.write(bs.read_text())
    bs.unlink()

    html = app.material('1MoS2-1')
    assert '1.24' in html  # Bader charge

    (tmp_path / 'AB2/1MoS2/1/bader.json').unlink()
    (tmp_path / 'AB2/1MoS2/1/results-asr.shift.json').unlink()
    (tmp_path / 'AB2/1MoS2/1/results-asr.convex_hull.json').unlink()
    html = app.material('1MoS2-1')
    assert '1.24' not in html  # Bader charge
