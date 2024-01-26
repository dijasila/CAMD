from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator

from camdweb.cli import main


def test_cli(tmp_path):
    f = tmp_path / 'H2'
    f.mkdir()
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    h2 = f / 'h2.xyz'
    atoms.write(h2)

    f = tmp_path / 'H'
    f.mkdir()
    atoms = Atoms('H')
    atoms.center(vacuum=2)
    atoms.calc = SinglePointCalculator(atoms, magmom=1.0)
    h = f / 'h.xyz'
    atoms.write(h)

    f = tmp_path / 'Cu'
    f.mkdir()
    b = 3.6 / 2
    atoms = Atoms('Cu', cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
    atoms.calc = EMT()
    atoms.get_forces()
    atoms.get_stress()
    cu = f / 'cu.xyz'
    atoms.write(cu)

    app = main([str(x) for x in [h2, h, cu]], run=False)
    assert 'Volume' in app.index()
