import os
from pathlib import Path

import pytest
from boddle import boddle
from camdweb.cmr.app import main
from camdweb.cmr.projects import abs3_bs, create_project_description, projects
from camdweb.test.cmr import create_db_file


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

    app.overview()
    app.favicon()

    papp = app.project_apps[name]

    papp.index_page()
    for material in papp.materials:
        papp.material_page(material.uid)

    if name == 'oqmd123':
        xyz = papp.download('id-1', 'xyz')
        assert 'energy=-27.0' in xyz

    if name == 'abs3':
        papp.png('abs3/bs-1.png')

    # Test also when png-files have already been generated:
    if name == 'abs3':
        papp.material_page('1')
    if name == 'lowdim':
        papp.material_page('a1')

    if name == 'ads1d':
        with boddle(query={'name': 'atoms', 'uid': 'id-1', 'data': '1'}):
            dct = papp.callback()
        assert 'data' in dct

    if name == 'abx2':
        mat = papp.materials['1']
        gap = mat.KS_gap
        assert isinstance(gap, float)
        with pytest.raises(AttributeError):
            mat.e_hull


def test_pd():
    create_project_description('hello')


def test_abs3():
    assert not abs3_bs({}, Path())
    assert not abs3_bs({'X': [1, 2], 'names': 'ABC'}, Path())


def test_pickle(in_tmp_path, tmp_path):
    create_db_file('abs3', tmp_path)
    app1 = main(['abs3.db', '--pickle'])
    (tmp_path / 'abs3.db').unlink()
    assert (tmp_path / 'abs3.pckl').is_file()
    app2 = main(['abs3.db'])
    n1 = len(app1.project_apps['abs3'].materials)
    n2 = len(app2.project_apps['abs3'].materials)
    assert n1 == n2 == 2
