from camdweb.materials import Material
from camdweb.panels.emass import EmassPanel
import json
from ase import Atoms
import numpy as np


def test_emass(tmp_path):
    atoms = Atoms()
    material = Material('x', atoms, folder=tmp_path)
    vbm_dct = {}
    vbm_dct['warping'] = 0
    vbm_dct['coords'] = [0.5, 0.0]
    vbm_dct['min_emass'] = 1.0
    vbm_dct['max_emass'] = 2.0
    vbm_dct['min_emass_direction'] = [1, 0]
    vbm_dct['max_emass_direction'] = [0, 1]
    vbm_dct['m_dos'] = np.sqrt(2)
    x = np.linspace(-1, 1, 10)
    X, Y = np.meshgrid(x, x)
    Z = np.sin(X * Y)
    vbm_dct['X'] = X.tolist()
    vbm_dct['Y'] = Y.tolist()
    vbm_dct['Z'] = Z.tolist()
    vbm_dct['barrier_found'] = False
    vbm_dct['dist_to_barrier'] = 10
    vbm_dct['extremum_depth'] = 10

    cbm_dct = vbm_dct

    dct = {}
    dct['vbm'] = vbm_dct
    dct['cbm'] = cbm_dct
    (tmp_path / 'emass.json').write_text(json.dumps(dct))

    emasspanel = EmassPanel()
    emasspanel.get_data(material)

    return
