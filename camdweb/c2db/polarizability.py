"""Optical polarizability."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from camdweb.html import image, table
from camdweb.panels.panel import Panel, PanelData, SkipPanel

HTML = """
<div class="row">
  <div class="col-6">
    <div id='bandstructure'></div>
  </div>
  <div class="col-6">
    {dos}
    {bz}
  </div>
</div>
"""
HTML = """
<div class="row">
  <div class="col-6">
   {x}
   {z}
  </div>
  <div class="col-6">
   {y}
   {table}
  </div>
</div>
"""

INFO = """\
The frequency-dependent polarisability in the long wave length limit (q=0)
calculated in the random phase approximation (RPA) without spin–orbit
interactions.  For metals a Drude term accounts for intraband transitions.
The contribution from polar lattice vibrations is added (see infrared
polarisability) and may be visible at low frequencies.
"""


class OpticalPolarizability(Panel):
    info = INFO

    def get_data(self, material):
        file = material.folder / 'polarizability.json'
        if not file.is_file():
            raise SkipPanel

        png = material.folder / 'rpa-pol-x.png'
        if not png.is_file():
            data = json.loads(file.read_text())
            omega_w, alpha_vw = json_dct_to_alpha_vw(data)
            create_optical_figures(
                omega_w,
                alpha_vw,
                material.plasmafrequency_x,
                material.plasmafrequency_y,
                material.folder)

        tab = table(
            ['Properties', ' '],
            self.table_rows(material,
                            [f'alpha{v}_el' for v in 'xyz'] +
                            [f'plasmafrequency_{v}]' for v in 'xy']))
        x, y, z = (image(material.folder / f'rpa-pol-{v}.png') for v in 'xyz')
        return PanelData(
            HTML.format(x=x, y=y, z=z, table=tab),
            title='Optical polarizability')


def json_dct_to_alpha_vw(data: dict) -> tuple(np.ndarray, np.ndarray):
    return (np.array(data['omega_w']),
            np.array([np.array(data[f'alpha{v}_re_w']) +
                      1j * np.array(data[f'alpha{v}_im_w'])
                      for v in 'xyz']))


def create_optical_figures(frequencies: np.ndarray,
                           alpha_vw: np.ndarray,
                           plasmafrequency_x: float,
                           plasmafrequency_y: float,
                           folder: Path) -> None:
    i2 = abs(frequencies - 50.0).argmin()
    frequencies = frequencies[:i2]
    alphax_w, alphay_w, alphaz_w = alpha_vw[:, :i2]

    ax = plt.figure().add_subplot(111)
    try:
        wpx = plasmafrequency_x
        if wpx > 0.01:
            alphaxfull_w = (alphax_w -
                            wpx**2 / (2 * np.pi * (frequencies + 1e-9)**2))
            ax.plot(
                frequencies,
                np.real(alphaxfull_w),
                '-',
                c='C1',
                label='real')
            ax.plot(
                frequencies,
                np.real(alphax_w),
                '--',
                c='C1',
                label='real (interband)')
        else:
            ax.plot(frequencies, np.real(alphax_w), c='C1', label='real')
    except AttributeError:
        ax.plot(frequencies, np.real(alphax_w), c='C1', label='real')
    ax.plot(frequencies, np.imag(alphax_w), c='C0', label='imag')

    fx, fy, fz = (folder / f'rpa-pol-{v}.png' for v in 'xyz')

    plot_polarizability(ax, frequencies, alphax_w, filename=fx, direction='x')

    ax = plt.figure().add_subplot(111)
    try:
        wpy = plasmafrequency_y
        if wpy > 0.01:
            alphayfull_w = (alphay_w -
                            wpy**2 / (2 * np.pi * (frequencies + 1e-9)**2))
            ax.plot(
                frequencies,
                np.real(alphayfull_w),
                '-',
                c='C1',
                label='real')
            ax.plot(
                frequencies,
                np.real(alphay_w),
                '--',
                c='C1',
                label='real (interband)')
        else:
            ax.plot(frequencies, np.real(alphay_w), c='C1', label='real')
    except AttributeError:
        ax.plot(frequencies, np.real(alphay_w), c='C1', label='real')

    ax.plot(frequencies, np.imag(alphay_w), c='C0', label='imag')
    plot_polarizability(ax, frequencies, alphay_w, filename=fy, direction='y')

    ax3 = plt.figure().add_subplot(111)
    ax3.plot(frequencies, np.real(alphaz_w), c='C1', label='real')
    ax3.plot(frequencies, np.imag(alphaz_w), c='C0', label='imag')
    plot_polarizability(ax3, frequencies, alphaz_w, filename=fz, direction='z')


def ylims(ws, data, wstart=0.0):
    i = abs(ws - wstart).argmin()
    x = data[i:]
    x1, x2 = x.real, x.imag
    y1 = min(x1.min(), x2.min()) * 1.02
    y2 = max(x1.max(), x2.max()) * 1.02
    return y1, y2


def plot_polarizability(ax, frequencies, alpha_w, filename, direction):
    import matplotlib.pyplot as plt
    ax.set_title(f'Polarization: {direction}')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel(r'Polarizability [$\mathrm{\AA}$]')
    ax.set_ylim(ylims(ws=frequencies, data=alpha_w, wstart=0.5))
    ax.legend()
    ax.set_xlim((0, 10))
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(filename)
    plt.close()


class IRPolarizability(Panel):
    def get_data(self, material):
        file = material.folder / 'ir-polarizability.json'
        if not file.is_file():
            raise SkipPanel

        png = material.folder / 'ir-pol-x.png'
        if not png.is_file():
            opt_data = json.loads(
                (material.folder / 'opt-polarizability.json').read_text())
            omega_el_w, alpha_el_vw = json_dct_to_alpha_vw(opt_data)
            data = json.loads(file.read_text())
            create_ir_figures(
                np.array(data['omega_w']),
                np.array(data['alpha_wvv']),
                omega_el_w,
                alpha_el_vw,
                data['maxphononfreq'],
                material.folder)

        tab = table(
            ['Properties', ' '],
            self.table_rows(material,
                            [f'alpha{v}_el' for v in 'xyz'] +
                            [f'plasmafrequency_{v}]' for v in 'xy']))

        x, y, z = (image(material.folder / f'ir-pol-{v}.png') for v in 'xyz')

        return PanelData(
            HTML.format(x=x, y=y, z=z, table=tab),
            title='Infrared polarizability')


def create_ir_figures(omega_w: np.ndarray,
                      alpha_wvv: np.ndarray,
                      omega_el_w: np.ndarray,
                      alpha_el_vw: np.ndarray,
                      maxphononfreq: float,
                      folder: Path) -> None:

    for v, (axisname, alpha_el_w) in enumerate(zip('xyz', alpha_el_vw)):
        create_plot_simple(
            ndim=2,
            maxomega=maxphononfreq * 1.5 * 1e3,
            omega_w=omega_w * 1e3,
            alpha_w=alpha_el_w,
            alphavv_w=alpha_wvv[:, v, v],
            omegatmp_w=omega_el_w * 1e3,
            axisname=axisname,
            fname=folder / f'ir-pol-{axisname}.png')


def create_plot_simple(*,
                       ndim,
                       omega_w, fname, maxomega, alpha_w,
                       alphavv_w, axisname,
                       omegatmp_w):
    from scipy.interpolate import interp1d

    re_alpha = interp1d(omegatmp_w, alpha_w.real)
    im_alpha = interp1d(omegatmp_w, alpha_w.imag)
    a_w = (re_alpha(omega_w) + 1j * im_alpha(omega_w) + alphavv_w)

    if ndim == 3:
        ylabel = r'Dielectric function'
        yvalues = 1 + 4 * np.pi * a_w
    else:
        power_txt = {2: '', 1: '^2', 0: '^3'}[ndim]
        unit = rf"$\mathrm{{\AA}}{power_txt}$"
        ylabel = rf'Polarizability [{unit}]'
        yvalues = a_w

    return mkplot(yvalues, axisname, fname, maxomega, omega_w, ylabel)


def mkplot(a_w, axisname, fname, maxomega, omega_w, ylabel):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(omega_w, a_w.real, c='C1', label='real')
    ax.plot(omega_w, a_w.imag, c='C0', label='imag')
    ax.set_title(f'Polarization: {axisname}')
    ax.set_xlabel('Energy [meV]')
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, maxomega)
    ax.legend()
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    return fname
