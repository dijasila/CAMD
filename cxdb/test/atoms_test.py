from ase import Atoms

from cxdb.atoms import AtomsPanel, plot_atoms


def test_1d():
    assert AtomsPanel(1).column_names['length'] == 'Length [Å]'


def test_plot():
    plot_atoms(Atoms('H', [[1, 1, 1]], cell=[2, 2, 2]))
