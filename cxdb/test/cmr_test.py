import os

from cxdb.cmr.app import main
from cxdb.test.cmr import create_db_files
from cxdb.cmr.projects import create_project_description


def test_cmr(tmp_path):
    os.chdir(tmp_path)
    project_descriptions = create_db_files(tmp_path)

    app = main([f'{name}.db' for name in project_descriptions])

    app.overview()

    app.index1('solar')
    app.material('solar', '1')

    html = app.material('oqmd123', 'id-1')
    assert 'http://oqmd.org' in html

    app.download_db_file('abs3')

    app.favicon()

    dct = app.callback('ads1d', {'name': 'atoms', 'uid': 'id-1', 'data': '1'})
    assert 'data' in dct

    gap = app.project_apps['abx2'].materials['1'].KS_gap
    assert isinstance(gap, float)


def test_pd():
    create_project_description('hello')
