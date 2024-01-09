from pathlib import Path

from ase import Atoms
from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.session import Session


def test_mat():
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    materials = Materials(
        [Material(Path(), 'x', atoms)],
        [AtomsPanel(3)])
    s = Session(1, ['uid'])
    rows, header, pages, new_columns = materials.get_rows(s)
    assert len(rows) == 1
    s.update('volume>1,stoichiometry=A', {})
    rows, header, pages, new_columns = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=A', {})
    rows, header, pages, new_columns = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=AB', {})
    rows, header, pages, new_columns = materials.get_rows(s)
    assert len(rows) == 0
