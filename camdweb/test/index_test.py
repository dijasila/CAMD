import pytest
from camdweb.filter import Index
from camdweb.parse import parse


@pytest.fixture(scope='module')
def index():
    return Index(
        [({'H': 2}, {'gap': 10.0, 'i1': 42, 'b1': False, 's1': 'OK'}),
         ({'Si': 2}, {'gap': 1.1}),
         ({'Si': 4}, {'gap': 1.1})])


@pytest.mark.parametrize(
    'filter, result',
    [('gap > 5.0', {0}),
     ('gap <= 5.0', {1, 2}),
     ('gap != 5.0', {0, 1, 2}),
     ('gap = 1.1', {1, 2}),
     ('Si > 2', {2}),
     ('Si', {1, 2}),
     ('SiH2', set()),
     ('gap < 2.0', {1, 2}),
     ('Si, gap < 2.0', {1, 2}),
     ('s1=OK', {0}),
     ('s1=OK2', set()),
     ('s1!=OK', set()),
     ('s1!=OK2', {0}),
     ('i1=42', {0}),
     ('i1!=42', set()),
     ('i1<42', set()),
     ('i1<=42', {0}),
     ('i1<=43', {0}),
     ('i1>=43', set()),
     ('i1>=42', {0}),
     ('b1=False', {0}),
     ('b1=0', {0}),
     ('Source=COD', SyntaxError),
     ('b1235=0', set())] +
    [(f, set()) for f in ['N', 'N=1', 'N!=0', 'N>0', 'N>=1', 'N<0']] +
    [(f, {0, 1, 2}) for f in ['N=0', 'N!=1', 'N<1', 'N<=0', 'N<=1']])
def test_index(index: Index, filter: str, result: set[int]):
    if isinstance(result, set):
        func = parse(filter)
        assert func(index) == result
    else:
        with pytest.raises(result):
            func = parse(filter)


@pytest.fixture(scope='module')
def mos_index():
    return Index(
        [({'Mo': 1, 'S': 2}, {}),
         ({'Mo': 2, 'S': 4}, {}),
         ({'Mo': 4, 'S': 8}, {}),
         ({'Mo': 1, 'S': 1}, {})])


@pytest.mark.parametrize(
    'filter, result',
    [('MoS2', {0, 1, 2}),
     ('MoS', {3}),
     ('MoS3', set()),
     ('Mo2S4', {1, 2}),
     ('Mo', {0, 1, 2, 3})])
def test_mos_index(mos_index: Index, filter: str, result: set[int]):
    func = parse(filter)
    assert func(mos_index) == result


def test_parse_errors():
    with pytest.raises(SyntaxError):
        parse('(Si=2')
    with pytest.raises(SyntaxError):
        parse('Si=2;a=b')


def test_index2():
    with pytest.raises(TypeError):
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
