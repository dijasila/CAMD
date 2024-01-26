import pytest

from ase import Atoms
from ase.calculators.emt import EMT

from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.bader import BaderPanel
from camdweb.panels.dos import DOSPanel
from camdweb.material import Material, Materials
from camdweb.web import CAMDApp
from boddle import boddle


@pytest.fixture(scope='module')
def c2db(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp('tmp-c2db')
    path = tmp_path / 'H2'
    path.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.write(path / 'structure.xyz')
    (path / 'dos.png').write_text('DOS')
    (path / 'bader.json').write_text('{"charges": [1.23, 0.0]}')
    c2db = CAMDApp(
        Materials([Material.from_file(path / 'structure.xyz', 'h2')],
                  [AtomsPanel(), DOSPanel(), BaderPanel()]),
        {'uid', 'volume', 'formula'},
        tmp_path)
    return c2db


def test_query_h2(c2db):
    with boddle(query={'filter': 'H=2'}):
        out = c2db.index()
    assert 'H<sub>2' in out


def test_query_sid_h2(c2db):
    with boddle(query={'sid': '0', 'filter': 'H=3,energy=42.0'}):
        out = c2db.index()
    assert 'H<sub>2' not in out


def test_query_stoichiometry_h2(c2db):
    with boddle(query={'stoichiometry': 'A', 'nspecies': '1'}):
        out = c2db.index()
    assert 'H<sub>2' in out


def test_various_queries(c2db):
    for query in [{'toggle': 'volume'},
                  {'toggle': 'volume'},
                  {'sort': 'volume'},
                  {'sort': 'volume'},
                  {'page': '0'}]:
        with boddle(query=query):
            c2db.index()


def test_material(c2db):
    out = c2db.material('h2')
    assert 'Atoms' in out
    assert 'Density of states' in out
    assert '1.23' in out


def test_callback(c2db):
    with boddle(query={'name': 'atoms', 'uid': 'h2', 'data': '2'}):
        dct = c2db.callback()
    assert 'data' in dct


def test_png(c2db):
    c2db.png('H2/dos.png')


def test_help(c2db):
    c2db.help()


@pytest.mark.parametrize('fmt, ref_substring', [
    ('cif', b'data_image0'),
    ('json', '"energy": 1.4'),
    ('xyz', 'energy=1.4'),
])
def test_download(c2db, fmt, ref_substring):
    # In principle we should test that the downloaded file is the same
    # atoms object as the one we downloaded.
    data = c2db.download('h2', fmt)
    assert ref_substring in data
