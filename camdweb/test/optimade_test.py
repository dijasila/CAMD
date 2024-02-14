import pytest
from optimade.server.exceptions import BadRequest  # type: ignore
from camdweb.filter import Index

@pytest.fixture(scope='module')
def make_index_db(path):
    index = Index(
        [({'H': 2},
          {'energy' : -1.0, 'hmm': 'a123'}),
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
     ['a+b', BadRequest],
     ['~a', BadRequest]])
def test_select(index, query, result, parser):
    if isinstance(result, list):
        tree1 = parser.parse(query)
    else:
        with pytest.raises(result):
            tree1 = parser.parse(query)
        return
    tree = parse_lark_tree(tree1)
    print(tree)
    selection = indexdb.select(tree)
    rows = indexdb.execute(selection)
    assert rows == result
