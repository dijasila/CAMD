import pytest
from camdweb.filter import parse, Index


@pytest.fixture(scope='module')
def index():
    return Index(
        [('H', {'H': 2}, {'gap': 10.0, 'i1': 42, 'b1': False, 's1': 'OK'}),
         ('Si', {'Si': 2}, {'gap': 1.1}),
         ('Si', {'Si': 4}, {'gap': 1.1})])


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
     pytest.param('Source=COD', set(), marks=[pytest.mark.xfail])] +
    [(f, set()) for f in ['N', 'N=1', 'N!=0', 'N>0', 'N>=1', 'N<0']] +
    [(f, {0, 1, 2}) for f in ['N=0', 'N!=1', 'N<1', 'N<=0', 'N<=1']])
def test_index(index: Index, filter: str, result: set[int]):
    func = parse(filter)
    assert func(index) == result


@pytest.fixture(scope='module')
def mos_index():
    return Index(
        [('MoS2', {'Mo': 1, 'S': 2}, {}),
         ('MoS2', {'Mo': 2, 'S': 4}, {}),
         ('MoS2', {'Mo': 4, 'S': 8}, {}),
         ('MoS', {'Mo': 1, 'S': 1}, {})])


@pytest.mark.parametrize(
    'filter, result',
    [('MoS2', {0, 1, 2}),
     ('MoS', {3}),
     ('MoS3', set()),
     ('Mo2S4', {1, 2})])
def test_mos_index(mos_index: Index, filter: str, result: set[int]):
    func = parse(filter)
    assert func(mos_index) == result


def test_parse_errors():
    with pytest.raises(SyntaxError):
        parse('(Si=2')
    with pytest.raises(SyntaxError):
        parse('Si=2;a=b')


def test_index2():
    with pytest.raises(ValueError):
        Index([('H', {'H': 2}, {'gap': 1j})])
    i = Index([('H', {'H': 2}, {'gap': 'high'})])
    with pytest.raises(ValueError):
        i.key('gap', '<', 'high')
    assert i.key('xx', '=', 117) == set()

    i = Index([('H', {'H': 2}, {'gap': 2.2, 'na': 27})])
    with pytest.raises(SyntaxError):
        i.float_key('gap', '~', 27)
    with pytest.raises(AssertionError):
        i.integer_key('na', '~', 27)
