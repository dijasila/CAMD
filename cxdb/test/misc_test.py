from cxdb.filter import str2obj, bisect


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_bisect():
    assert bisect([0, 1, 2], -0.5) == 0
    assert bisect(range(100), 7.5) == 8
