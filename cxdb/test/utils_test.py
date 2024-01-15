from cxdb.utils import Range, RangeX


def test_range():
    r = Range('ABC', 'x')
    ','.join(r.get_filter_strings({})) == ''
    ','.join(r.get_filter_strings({'to_x': '1'})) == 'x<=1'
    ','.join(r.get_filter_strings({'from_x': '1'})) == 'x>=1'
    ','.join(r.get_filter_strings({'from_x': '1', 'to_x': '2'})) == 'x>=1,x<=2'


def test_rangeX():
    r = RangeX('ABC', 'x', ['A', 'B'])
    assert r.get_filter_strings(
        {'x': 'y', 'from_x': '0.5'}) == ['y>=0.5']
    assert r.get_filter_strings(
        {'x': 'y', 'to_x': '0.5'}) == ['y<=0.5']
