from cxdb.atoms import AtomsPanel


def test_1d():
    assert AtomsPanel(1).column_names['length'] == 'Length [Å]'
