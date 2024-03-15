import numpy as np
from asr.effective_masses import is_band_warped, map_to_IBZ
from asr.effective_masses import unit_to_electron_mass


def get_emass_data(data, atoms) -> dict:
    webpanel_data = {}
    band_names = ['vbm', 'cbm']
    if data:
        unit_cell = data['unit_cell']
        webpanel_data['unit_cell'] = unit_cell
        for band_name in band_names:
            if band_name + '_data' not in data:
                continue
            band_data = data[band_name + '_data']
            webpanel_band_data = {}
            for key in band_data:
                band_data[key] = np.asarray(band_data[key])

            band_warped = is_band_warped(band_data['iems_warping'])
            webpanel_band_data['warping'] = band_data['iems_warping']
            coords_ibz = map_to_IBZ(band_data['coords_cartesian'], atoms)
            coords_kbasis = coords_ibz @ unit_cell.T
            webpanel_band_data['coords'] = coords_kbasis.flatten()

            phi = band_data['iems_phi']
            iems_coefficients_k = band_data['iems_coefficients_ks'][:, 0]
            if band_warped:
                iems_coefficients_k = band_data['iems_coefficients_ks'][:, 0]
                inverse_max_emass = np.min(np.abs(iems_coefficients_k))
                inverse_min_emass = np.max(np.abs(iems_coefficients_k))
                if len(np.unique(np.sign(iems_coefficients_k))) > 1:
                    inverse_max_emass = 0

                max_emass_angle = phi[np.argmin(np.abs(iems_coefficients_k))]
                max_emass_direction = np.array(
                    [np.cos(max_emass_angle), np.sin(max_emass_angle)])
                min_emass_angle\
                    = phi[np.argmax(np.abs(iems_coefficients_k))]
                min_emass_direction = np.array(
                    [np.cos(min_emass_angle), np.sin(min_emass_angle)])
            else:
                eigvals = band_data['fit_eigvals']
                eigvecs = band_data['fit_eigvecs']
                max_emass_idx = np.argmin(abs(eigvals))
                min_emass_idx = (1 - max_emass_idx) % 2
                inverse_max_emass = np.abs(eigvals[max_emass_idx]) / 2
                inverse_min_emass = np.abs(eigvals[min_emass_idx]) / 2
                if np.sign(eigvals[0]) != np.sign(eigvals[1]):
                    inverse_max_emass = 0
                max_emass_direction = eigvecs[:, max_emass_idx]
                min_emass_direction = eigvecs[:, min_emass_idx]

            if inverse_min_emass > 0:
                min_emass = (1 / inverse_min_emass) * unit_to_electron_mass
            else:
                min_emass = np.inf

            if inverse_max_emass > 0:
                max_emass = (1 / inverse_max_emass) * unit_to_electron_mass
            else:
                max_emass = np.inf

            webpanel_band_data['min_emass'] = min_emass
            webpanel_band_data['max_emass'] = max_emass
            webpanel_band_data['min_emass_direction'] = min_emass_direction
            webpanel_band_data['max_emass_direction'] = max_emass_direction
            if band_warped:
                m_dos = band_data['iems_m_dos'] * unit_to_electron_mass
            else:
                if inverse_max_emass * inverse_min_emass > 0:
                    m_dos = unit_to_electron_mass\
                        / np.sqrt(inverse_max_emass * inverse_min_emass)
                else:
                    m_dos = np.inf
            webpanel_band_data['m_dos'] = m_dos
            X = band_data['contour_kx']
            Y = band_data['contour_ky']
            f0 = band_data['fit_f0']
            Z = band_data['contour_energies'] - f0

            webpanel_band_data['X'] = X
            webpanel_band_data['Y'] = Y
            webpanel_band_data['Z'] = Z
            energy_levels = band_data['barrier_levels']
            dx = X[0, 1] - X[0, 0]
            barrier_R = band_data['barrier_R']
            barrier_found = False
            if np.any(np.diff(barrier_R) > 2.9 * dx):
                discont_idx = np.nonzero(
                    np.diff(barrier_R) > 2.9 * dx)[0][0]
                dist_to_barrier = barrier_R[discont_idx]  # size of cbm in 1/Å
                barrier_found = True
            # depth of cbm in meV
                extremum_depth = energy_levels[discont_idx] * 1000

            else:
                dist_to_barrier = barrier_R.max()
                extremum_depth = energy_levels.max() * 1000

            webpanel_band_data['barrier_found'] = barrier_found
            webpanel_band_data['dist_to_barrier'] = dist_to_barrier
            webpanel_band_data['extremum_depth'] = extremum_depth

            webpanel_data[band_name] = webpanel_band_data
    return webpanel_data
