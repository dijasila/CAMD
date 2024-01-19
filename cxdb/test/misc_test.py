import pytest
from ase import Atoms
from cxdb.panels.asr_panel import thing2html
from cxdb.panels.bader import BaderPanel
from cxdb.filter import str2obj
from cxdb.material import Material


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_no_bader(tmp_path):
    material = Material(tmp_path, 'x1', Atoms())
    assert BaderPanel().get_html(material, None) == ('', '')


def test_thing():
    with pytest.raises(ValueError):
        thing2html({'type': '?'}, 'x')
