import pytest
from camdweb.filter import Index
from camdweb.optimade.filter import create_parse_function, select


@pytest.fixture(scope='module')
def parse():
    return create_parse_function()


@pytest.fixture(scope='module')
def index():
    index = Index(
        [({'H': 2},
          {'energy': -1.0, 'hmm': 'a123',
           'nspecies': 1, 'stoichiometry': 'A', 'formula': 'H2'}),
         ({'C': 1, 'O': 1},
          {'x': 7, 'hmm': 'hi!',
           'nspecies': 2, 'stoichiometry': 'AB', 'formula': 'CO'})])
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
     ['hmm="abc"', []],
     ['hmm="hi!" OR hmm="a123" AND energy=-1.0', [1, 2]],
     ['(hmm="hi!" OR hmm="a123") AND energy=-1.0', [1]],
     ['elements HAS "C"', [2]],
     ['elements HAS ALL "C"', [2]],
     ['elements HAS ANY "C"', [2]],
     ['NOT (elements HAS "C")', [1]],
     ['elements HAS "C" AND x=7', [2]],
     ['elements HAS "Cu"', []],
     ['a+b', ValueError],
     ['~a', ValueError]])
def test_select(index, query, result, parse):
    from lark.exceptions import UnexpectedCharacters
    if isinstance(result, list):
        tree = parse(query)
    else:
        with pytest.raises(UnexpectedCharacters):
            parse(query)
        return
    print(tree)
    selection = select(tree, index)
    print(selection)
    rows = sorted(i + 1 for i in selection)
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
