from ase import Atoms
from ase.calculators.emt import EMT

from cxdb.atoms import AtomsPanel
from cxdb.bader import BaderPanel
from cxdb.dos import DOSPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp


def test_app(tmp_path):
    f = tmp_path / 'H2'
    f.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
    atoms.center(vacuum=2)
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.write(f / 'structure.xyz')
    (f / 'dos.png').write_text('DOS')
    (f / 'bader.json').write_text('{"charges": [1.23, 0.0]}')
    c2db = CXDBApp(Materials([Material.from_file(f / 'structure.xyz', 'h2')],
                             [AtomsPanel(3), DOSPanel(), BaderPanel()]),
                   {'uid', 'volume', 'formula'},
                   tmp_path)
    out = c2db.index({'filter': 'H=2'})
    assert 'H<sub>2' in out
    out = c2db.index({'filter': 'H=3,energy=42.0'})
    assert 'H<sub>2' not in out
    out = c2db.index({'stoichiometry': 'A', 'nspecies': '1'})
    assert 'H<sub>2' in out
    c2db.index({'toggle': 'volume'})
    c2db.index({'toggle': 'volume'})
    c2db.index({'sort': 'volume'})
    c2db.index({'sort': 'volume'})
    c2db.index({'page': '0'})
    out = c2db.material('h2')
    assert 'Atoms' in out
    assert 'Density of states' in out
    assert '1.23' in out
    dct = c2db.callback(dict(name='atoms', uid='h2', data='2'))
    assert 'data' in dct
    c2db.png('h2', 'dos.png')
    c2db.help()
