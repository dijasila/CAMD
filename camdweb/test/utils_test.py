import pytest

from camdweb.html import Range, RangeS, RangeX
from camdweb.utils import process_pool


@pytest.fixture(scope='module')
def r():
    return Range('ABC', 'x')


@pytest.mark.parametrize(
    'query, result',
    [({}, ''),
     ({'to_x': '1'}, 'x<=1'),
     ({'from_x': '1'}, 'x>=1'),
     ({'from_x': '1', 'to_x': '2'}, 'x>=1,x<=2')])
def test_range(r, query, result):
    assert ','.join(r.get_filter_strings(query)) == result


def test_range_positive():
    r = Range('ABC', 'x', nonnegative=True)
    assert r.get_filter_strings({'from_x': '0'}) == []


def test_range_x():
    r = RangeX('ABC', 'x', ['A', 'B'])
    assert r.get_filter_strings(
        {'x': 'y', 'from_x': '0.5'}) == ['y>=0.5']
    assert r.get_filter_strings(
        {'x': 'y', 'to_x': '0.5'}) == ['y<=0.5']


def test_range_s():
    r = RangeS('ABC', 'x', ['A', 'B'])
    assert r.get_filter_strings({'from_x': '0.5'}) == ['x>=0.5']
    assert r.get_filter_strings({'to_x': '0.5'}) == ['x<=0.5']


def f(x):
    return x


@pytest.mark.parametrize('n', [1, 2])
def test_ppool(n):
    s = ''
    with process_pool(n) as pool:
        for x in pool.imap_unordered(f, 'abc'):
            s += x
    assert ''.join(sorted(s)) == 'abc'
