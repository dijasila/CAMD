from camdweb.materials import Material
from camdweb.panels.emass import EmassPanel
from camdweb.c2db.emass import get_emass_data
import json
from ase.build import bulk
from ase import Atoms
import numpy as np


def test_emass_data(tmp_path):
    atoms = bulk('Si')

    band_data = {}
    band_data['unit_cell'] = atoms.cell[:2, :2]
    vbm_data = {}
    vbm_data['fit_eigvals'] = np.array([0.1, 0.5])
    vbm_data['fit_eigvecs'] = np.array([[1, 0], [0, 1]])
    vbm_data['fit_f0'] = 0
    vbm_data['iems_phi'] = np.linspace(0, 2 * np.pi, 11)
    vbm_data['iems_coefficients_ks'] = np.ones((11, 3))
    vbm_data['coords_cartesian'] = np.array([0, 0.1])
    vbm_data['iems_warping'] = 0
    vbm_data['iems_m_dos'] = np.sqrt(0.1 * 0.5)
    x = np.linspace(-1, 1, 10)
    X, Y = np.meshgrid(x, x)
    Z = np.sin(X * Y)
    vbm_data['contour_kx'] = X
    vbm_data['contour_ky'] = Y
    vbm_data['contour_energies'] = Z
    vbm_data['barrier_levels'] = np.array([0, 0.1, 0.2])
    vbm_data['barrier_R'] = np.array([0, 0.01, 0.02])

    cbm_data = vbm_data.copy()
    cbm_data['iems_warping'] = 1.0

    band_data['vbm_data'] = vbm_data
    band_data['cbm_data'] = cbm_data

    webpanel_data = get_emass_data(band_data, atoms)

    webpanel_vbmdata = webpanel_data['vbm']
    keys = {'warping', 'coords', 'min_emass', 'max_emass',
            'min_emass_direction', 'max_emass_direction', 'm_dos',
            'X', 'Y', 'Z', 'barrier_found', 'dist_to_barrier',
            'extremum_depth'}
    for key in webpanel_vbmdata:
        assert key in keys

    return


def test_emass_panel(tmp_path):
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
