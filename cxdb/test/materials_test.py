from pathlib import Path

from ase import Atoms
from cxdb.panels.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.session import Session


def test_mat():
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    materials = Materials(
        [Material(Path(), 'x', atoms)],
        [AtomsPanel()])
    s = Session(1, ['uid'])
    rows, header, pages, new_columns = materials.get_rows(s)
    assert (rows, header, pages, new_columns) == (
        [('x', ['x'])],
        [('uid', 'Unique ID')],
        [(0, 'previous'), (0, 'next'), (0, '1-1')],
        [('formula', 'Formula'),
         ('stoichiometry', 'Stoichiometry'),
         ('nspecies', 'Number of species'),
         ('volume', 'Volume [Å<sup>3</sup>]')])
    s.update('volume>1,stoichiometry=A', {})
    rows, _, _, _ = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=A', {})
    rows, _, _, _ = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=AB', {})
    rows, _, _, _ = materials.get_rows(s)
    assert len(rows) == 0
