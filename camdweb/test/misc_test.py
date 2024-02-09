import pytest
from camdweb.c2db.asr_panel import thing2html
from camdweb.parse import str2obj
from camdweb.oqmd12345.app import oqmd
from camdweb.utils import cod, icsd


def test_str2obj():
    assert str2obj('True') is True
    assert str2obj('False') is False


def test_thing():
    with pytest.raises(ValueError):
        thing2html({'type': '?'}, 'x')


def test_oqmd():
    oqmd(27)


def test_formatters():
    cod(117, True)
    icsd('117', True)
