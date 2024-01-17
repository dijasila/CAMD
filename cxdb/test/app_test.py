from ase import Atoms
from ase.calculators.emt import EMT

from cxdb.panels.atoms import AtomsPanel
from cxdb.panels.bader import BaderPanel
from cxdb.panels.dos import DOSPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp
from boddle import boddle


def test_app(tmp_path):
    f = tmp_path / 'H2'
    f.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    atoms.calc = EMT()
    atoms.get_potential_energy()
    atoms.write(f / 'structure.xyz')
    (f / 'dos.png').write_text('DOS')
    (f / 'bader.json').write_text('{"charges": [1.23, 0.0]}')
    c2db = CXDBApp(Materials([Material.from_file(f / 'structure.xyz', 'h2')],
                             [AtomsPanel(), DOSPanel(), BaderPanel()]),
                   {'uid', 'volume', 'formula'},
                   tmp_path)
    with boddle(query={'filter': 'H=2'}):
        out = c2db.index()
    assert 'H<sub>2' in out

    with boddle(query={'sid': '0', 'filter': 'H=3,energy=42.0'}):
        out = c2db.index()
    assert 'H<sub>2' not in out

    with boddle(query={'stoichiometry': 'A', 'nspecies': '1'}):
        out = c2db.index()
    assert 'H<sub>2' in out

    for query in [{'toggle': 'volume'},
                  {'toggle': 'volume'},
                  {'sort': 'volume'},
                  {'sort': 'volume'},
                  {'page': '0'}]:
        with boddle(query=query):
            c2db.index()

    out = c2db.material('h2')
    assert 'Atoms' in out
    assert 'Density of states' in out
    assert '1.23' in out

    with boddle(query={'name': 'atoms', 'uid': 'h2', 'data': '2'}):
        dct = c2db.callback()
    assert 'data' in dct

    c2db.png('h2', 'dos.png')

    c2db.help()
