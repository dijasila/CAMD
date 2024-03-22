import numpy as np
import pytest
from ase.lattice import BCT

from camdweb.bandstructure import PlotUtil, prettify_labels


def test_plotutil():
    bct = BCT(3, 5)
    nkpts = 123
    bandpath = bct.bandpath(npoints=nkpts)
    bs = bandpath.free_electron_band_structure()
    nbands = bs.energies.shape[-1]

    assert bs.energies.shape == (1, nkpts, nbands)

    energies_1kn = bs.energies
    assert energies_1kn.shape[0] == 1
    energies_skn = bs.energies.repeat(2, axis=0)
    assert energies_skn.shape[0] == 2
    energies_skn[1] *= 1.02

    energy_soc_mk = energies_skn[0] * 1.04
    spin_projection_mk = np.sin(energy_soc_mk)

    util = PlotUtil(
        path=bandpath,
        energy_nosoc_skn=energies_skn,
        energy_soc_mk=energy_soc_mk,
        spin_zprojection_soc_mk=spin_projection_mk,
        fermilevel=2.0,
        emin=-1.0,
        emax=10.0,
        spin_axisname='ø')

    # Can we assert something?
    util.plot()


def test_kpoint_labels():
    labels = ['G', 'B', 'C', 'D', 'EFG123']
    coords = [1.0, 2.0, 3.0, 3.0, 4.0]  # discontinuous bandpath ABC|DE
    newlabels, newcoords = prettify_labels(labels, coords)
    assert newcoords == pytest.approx([1.0, 2.0, 3.0, 4.0])
    assert newlabels == ['Γ', 'B', 'C,D', 'EFG<sub>123</sub>']
