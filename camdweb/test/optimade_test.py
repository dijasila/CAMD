import pytest
from camdweb.filter import Index
from camdweb.optimade.filter import create_parse_function


@pytest.fixture(scope='module')
def parse():
    return create_parse_function()


@pytest.fixture(scope='module')
def index():
    index = Index(
        [({'H': 2},
          {'energy': -1.0, 'hmm': 'a123'}),
         ({'C': 1, 'O': 1},
          {'x': 7, 'hmm': 'hi!'})])
    return index


@pytest.mark.parametrize(
    'query, result',
    [['x=7', [2]],
     ['y=7', []],
     ['x=8', []],
     ['8=x', []],
     ['nelements<2', [1]],
     ['2>nelements', [1]],
     ['chemical_formula_reduced != "H"', [2]],
     ['"H" != chemical_formula_reduced', [2]],
     ['chemical_formula_anonymous != "A"', [2]],
     ['"A" != chemical_formula_anonymous', [2]],
     ['chemical_formula_reduced = "H"', [1]],
     ['"H" = chemical_formula_reduced', [1]],
     ['chemical_formula_anonymous = "A"', [1]],
     ['"A" = chemical_formula_anonymous', [1]],
     ['chemical_formula_hill CONTAINS "H"', [1]],
     ['chemical_formula_hill STARTS "H"', [1]],
     ['chemical_formula_hill ENDS "2"', [1]],
     ['id=1', [1]],
     ['hmm="abc"', []],
     ['hmm ENDS WITH "3"', [1]],
     ['hmm STARTS WITH "3"', []],
     ['hmm CONTAINS "i"', [2]],
     ['hmm="hi!" OR hmm="a123" AND energy=-1.0', [1, 2]],
     ['(hmm="hi!" OR hmm="a123") AND energy=-1.0', [1]],
     ['structure_features HAS "bulk"', []],
     ['NOT structure_features HAS "assemblies"', [1, 2]],
     ['elements HAS "C"', [2]],
     ['elements HAS ALL "C"', [2]],
     ['elements HAS ANY "C"', [2]],
     ['NOT (elements HAS "C")', [1]],
     ['elements HAS "C" AND x=7', [2]],
     ['elements HAS "Cu"', []],
     ['a+b', ValueError],
     ['~a', ValueError]])
def test_select(index, query, result, parse):
    if isinstance(result, list):
        tree = parse(query)
    else:
        with pytest.raises(result):
            parse(query)
        return
    print(tree)
    selection = index.select(tree)
    rows = index.execute(selection)
    assert rows == result


@pytest.mark.parametrize(
    'query, expected',
    [('x HAS ALL "Au", "Cu"',
      [('IDENTIFIER', 'x'),
       [('HAS', 'HAS'), ('ALL', 'ALL'), ['Au', 'Cu']]]),
     ('x=1 AND y=2',
      ('AND',
       [[('IDENTIFIER', 'x'), [('OPERATOR', '='), 1]],
        [('IDENTIFIER', 'y'), [('OPERATOR', '='), 2]]])),
     ('x=1 OR y=2',
      ('OR',
       [[('IDENTIFIER', 'x'), [('OPERATOR', '='), 1]],
        [('IDENTIFIER', 'y'), [('OPERATOR', '='), 2]]])),
     ('NOT x=1',
      [('NOT', 'NOT'),
       [('IDENTIFIER', 'x'), [('OPERATOR', '='), 1]]])])
def test_parser(query, expected, parse):
    tree = parse(query)
    print(tree)
    assert tree == expected


def test_create_lark_parser(parse):
    assert parse('x=2') == [('IDENTIFIER', 'x'), [('OPERATOR', '='), 2]]
