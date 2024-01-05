from cxdb.filter import str2obj, bisect
from cxdb.bader import BaderPanel
from cxdb.material import Material
from ase import Atoms


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_bisect():
    assert bisect([0, 1, 2], -0.5) == 0
    assert bisect(range(100), 7.5) == 8


def test_no_bader(tmp_path):
    material = Material(tmp_path, 'x1', Atoms())
    assert BaderPanel().get_html(material, None) == ('', '')
