import pytest
from camdweb.oqmd12345.app import main

DATA = """{
"oqmd_entry_id": 588771,
"etot": -0.014277391139242677,
"gap": 8.620173558602595,
"gap_dir": 8.620173558602595,
"magstate": "NM",
"fmax": 2.2102264799446573e-29,
"smax": 0.003093088350535923
}
"""

XYZ = ('1\n'
       'Lattice="3.8 0 0 1.9 3.3 0 1.9 1.1 3.1" '
       'Properties=species:S:1:pos:R:3:initial_magmoms:R:1 pbc="T T T"\n'
       'Ar       0.0       0.0       0.0       0.0\n')


@pytest.fixture
def oqmd12345(tmp_path_factory):
    root = tmp_path_factory.mktemp('tmp-oqmd')
    path = root / 'A/1Ar/588771'
    path.mkdir(parents=True)
    (path / 'data.json').write_text(DATA)
    (path / 'structure.xyz').write_text(XYZ)
    return root


def test_everything(oqmd12345):
    app = main(oqmd12345)
    assert len(app.materials) == 1
    app.index_page()
    html = app.material_page('1Ar-588771')
    assert '8.62' in html  # band gap
