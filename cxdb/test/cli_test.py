from ase import Atoms

from cxdb.cli import main


def test_cli(tmp_path):
    f = tmp_path / 'H2'
    f.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    h2 = f / 'h2.xyz'
    atoms.write(h2)
    app = main([str(h2)], run=False)
    assert 'Volume' in app.index()
