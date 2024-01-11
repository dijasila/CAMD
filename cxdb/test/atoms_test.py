from pathlib import Path

from ase import Atoms

from cxdb.material import Material
from cxdb.panels.atoms import AtomsPanel, plot_atoms


def test_1d():
    ap = AtomsPanel()
    mat = Material(Path(), 'x', Atoms('H', [[2.5, 2.5, 0]],
                                      cell=[5, 5, 1],
                                      pbc=[False, False, True]))
    ap.update_data(mat)
    assert ap.column_names['length'] == 'Length [Å]'


def test_plot():
    plot_atoms(Atoms('H', [[1, 1, 1]], cell=[2, 2, 2]))
