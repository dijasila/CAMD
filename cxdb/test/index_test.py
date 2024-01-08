import pytest
from cxdb.filter import parse, Index


def test_index():
    i = Index([({'H': 2}, {'gap': 10.0, 'i1': 42, 'b1': False, 's1': 'OK'}),
               ({'Si': 2}, {'gap': 1.1}),
               ({'Si': 4}, {'gap': 1.1})])
    f = parse('gap > 5.0')
    assert f(i) == {0}
    f = parse('gap <= 5.0')
    assert f(i) == {1, 2}
    f = parse('gap != 5.0')
    assert f(i) == {0, 1, 2}
    f = parse('gap = 1.1')
    assert f(i) == {1, 2}
    f = parse('Si > 2')
    assert f(i) == {2}
    f = parse('Si')
    assert f(i) == {1, 2}
    f = parse('SiH2')
    assert f(i) == set()
    f = parse('gap < 2.0')
    assert f(i) == {1, 2}
    f = parse('Si, gap < 2.0')
    assert f(i) == {1, 2}

    f = parse('s1=OK')
    assert f(i) == {0}
    f = parse('s1=OK2')
    assert f(i) == set()
    f = parse('s1!=OK')
    assert f(i) == set()
    f = parse('s1!=OK2')
    assert f(i) == {0}

    f = parse('i1=42')
    assert f(i) == {0}
    f = parse('i1!=42')
    assert f(i) == set()
    f = parse('i1<42')
    assert f(i) == set()
    f = parse('i1<=42')
    assert f(i) == {0}
    f = parse('i1<=43')
    assert f(i) == {0}

    f = parse('b1=False')
    assert f(i) == {0}
    f = parse('b1=0')
    assert f(i) == {0}

    for x in ['N', 'N=1', 'N!=0', 'N>0', 'N>=1', 'N<0']:
        f = parse(x)
        assert f(i) == set()
    for x in ['N=0', 'N!=1', 'N<1', 'N<=0', 'N<=1']:
        f = parse(x)
        assert f(i) == {0, 1, 2}

    with pytest.raises(SyntaxError):
        parse('(Si=2')
    with pytest.raises(SyntaxError):
        parse('Si=2;a=b')


def test_index2():
    with pytest.raises(ValueError):
        Index([({'H': 2}, {'gap': 1j})])
    i = Index([({'H': 2}, {'gap': 'high'})])
    with pytest.raises(ValueError):
        i.key('gap', '<', 'high')
    assert i.key('xx', '=', 117) == set()

    i = Index([({'H': 2}, {'gap': 2.2, 'na': 27})])
    with pytest.raises(SyntaxError):
        i.float_key('gap', '~', 27)
    with pytest.raises(AssertionError):
        i.integer_key('na', '~', 27)
