import pickle
from pathlib import Path

import pytest
from ase import Atoms
from camdweb.material import Material, Materials
from camdweb.c2db.asr_panel import Row
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.session import Session


@pytest.fixture(scope='module')
def material():
    atoms = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], pbc=True)
    atoms.center(vacuum=2)
    material = Material(Path(), 'x', atoms)
    return material


def test_materials(material):
    materials = Materials([material], [AtomsPanel()])
    s = Session(1, ['uid'])
    rows, header, pages, new_columns, error = materials.get_rows(s)
    print(new_columns)
    assert (rows, header, pages, new_columns, error) == (
        [('x', ['x'])],
        [('uid', 'Unique ID')],
        [(0, '«'), (0, '<'), (0, '1-1'), (0, '>'), (0, '»')],
        [('formula', 'Formula'),
         ('reduced', 'Reduced formula'),
         ('stoichiometry', 'Stoichiometry'),
         ('nspecies', 'Number of species'),
         ('natoms', 'Number of atoms'),
         ('volume', 'Unit cell volume [Å<sup>3</sup>]')],
        '')
    s.update('volume>1,stoichiometry=A', {})
    rows, _, _, _, _ = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=A', {})
    rows, _, _, _, _ = materials.get_rows(s)
    assert len(rows) == 1
    s.update('stoichiometry=AB', {})
    rows, _, _, _, _ = materials.get_rows(s)
    assert len(rows) == 0
    s.update('Ha>100', {})
    rows, _, _, _, error = materials.get_rows(s)
    assert len(rows) == 0
    assert error == 'Unknown chemical symbol "Ha"'


def test_attribute_error(material):
    with pytest.raises(AttributeError):
        material.asdf


def test_pickle(material):
    pickle.loads(pickle.dumps(material))


def ________test_collision():
    class MyPanel(Panel):
        column_names = {'formula': '...'}

        def get_html(self, material, materials):
            return

    with pytest.raises(ValueError):
        panel = MyPanel()
        panel.get_html(1, 2)
        Materials([], [panel])


def test_row(material):
    assert Row(material).toatoms().pbc.all()
