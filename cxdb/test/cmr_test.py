import os
from pathlib import Path

import pytest
from cxdb.cmr.app import main
from cxdb.cmr.projects import abs3_bs, create_project_description
from cxdb.test.cmr import create_db_files


def test_cmr(tmp_path):
    os.chdir(tmp_path)
    project_descriptions = create_db_files(tmp_path)

    app = main([f'{name}.db' for name in project_descriptions])

    app.overview()

    for name in project_descriptions:
        app.index1(name)
        for material in app.project_apps[name].materials:
            app.material(name, material.uid)

    html = app.material('oqmd123', 'id-1')
    assert 'http://oqmd.org' in html

    app.download_db_file('abs3')
    app.png('abs3', '1')
    app.material('abs3', '1')

    app.favicon()

    dct = app.callback('ads1d', {'name': 'atoms', 'uid': 'id-1', 'data': '1'})
    assert 'data' in dct

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
