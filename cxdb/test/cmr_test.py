import os
from pathlib import Path

import pytest
from cxdb.cmr.app import main
from cxdb.cmr.projects import abs3_bs, create_project_description, projects
from cxdb.test.cmr import create_db_file
from boddle import boddle


@pytest.fixture
def in_tmp_path(tmp_path):
    orig_cwd = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(orig_cwd)


@pytest.mark.parametrize('project_name', projects)
def test_cmr(in_tmp_path, tmp_path, project_name):
    name = project_name
    create_db_file(project_name, tmp_path)

    app = main([f'{name}.db'])

    app.index1(name)
    for material in app.project_apps[name].materials:
        app.material(name, material.uid)

    app.overview()
    app.favicon()

    if name == 'oqmd123':
        html = app.material('oqmd123', 'id-1')
        assert 'http://oqmd.org' in html

    if name == 'abs3':
        app.download_db_file('abs3')
        app.png('abs3', '1')
        # Test also when png-files have already been generated:
        app.material('abs3', '1')

    if name == 'lowdim':
        app.material('lowdim', 'a1')

    if name == 'ads1d':
        with boddle(query={'name': 'atoms', 'uid': 'id-1', 'data': '1'}):
            dct = app.callback('ads1d')
        assert 'data' in dct

    if name == 'abx2':
        mat = app.project_apps['abx2'].materials['1']
        gap = mat.KS_gap
        assert isinstance(gap, float)
        with pytest.raises(AttributeError):
            mat.E_hull


def test_pd():
    create_project_description('hello')


def test_abs3():
    assert not abs3_bs({}, Path())
    assert not abs3_bs({'X': [1, 2], 'names': 'ABC'}, Path())
