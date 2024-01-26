import pytest
from ase import Atoms
from camdweb.c2db.asr_panel import thing2html
from camdweb.panels.bader import BaderPanel
from camdweb.filter import str2obj
from camdweb.material import Material


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_no_bader(tmp_path):
    material = Material(tmp_path, 'x1', Atoms())
    assert BaderPanel().get_html(material, None) == ('', '')


def test_thing():
    with pytest.raises(ValueError):
        thing2html({'type': '?'}, 'x')
