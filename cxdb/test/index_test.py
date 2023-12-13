from cxdb.query import parse, Index


def test_index():
    i = Index([({'H': 2}, {'gap': 10.0}),
               ({'Si': 2}, {'gap': 1.1}),
               ({'Si': 4}, {'gap': 1.1})])
    f = parse('gap > 5.0')
    assert f(i) == {0}
    f = parse('Si > 2')
    assert f(i) == {2}
    f = parse('Si')
    assert f(i) == {1, 2}
    f = parse('gap < 2.0')
    assert f(i) == {1, 2}
    f = parse('Si, gap < 2.0')
    assert f(i) == {1, 2}
