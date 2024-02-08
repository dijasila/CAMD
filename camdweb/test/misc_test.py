import pytest
from camdweb.c2db.asr_panel import thing2html
from camdweb.panels.bader import BaderPanel
from camdweb.parse import str2obj
from camdweb.material import Material


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_no_bader():
    material = Material('x1')
    with pytest.raises(StopIteration):
        next(BaderPanel().get_html(material))


def test_thing():
    with pytest.raises(ValueError):
        thing2html({'type': '?'}, 'x')
