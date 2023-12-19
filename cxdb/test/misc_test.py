from cxdb.filter import str2obj


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False
