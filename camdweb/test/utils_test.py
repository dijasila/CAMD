from camdweb.html import Range, RangeX, RangeS
import pytest


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
