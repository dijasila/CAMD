from cxdb.web import C2DB
from cxdb.material import Material
from ase import Atoms
from ase.calculators.emt import EMT

h2 = """2

H 0.0 0.0 0.0
H 0.7 0.0 0.0
"""


def test_app(tmp_path):
    f = tmp_path / 'H2'
    f.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
    atoms.center(vacuum=2)
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.write(f / 'rlx.traj')
    (f / 'dos.png').write_text('DOS')
    (f / 'bader.json').write_text('{"charges": [1.23, 0.0]}')
    c2db = C2DB({'h2': Material(f, 'h2')})
    out = c2db.index()
    assert 'H<sub>2' in out
    out = c2db.material('h2')
    assert 'Atoms' in out
    assert 'Density of states' in out
    assert '1.23' in out
