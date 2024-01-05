import os

from cxdb.c2db import copy_materials, main
from cxdb.test.c2db import create_tree


def test_c2db(tmp_path):
    create_tree(tmp_path)
    os.chdir(tmp_path)
    copy_materials(tmp_path, ['MoS2*'])
    app = main(tmp_path)
    app.index()
    print(app.materials._materials)
    html = app.material('1MoS2-1')
    assert '15.000' in html
