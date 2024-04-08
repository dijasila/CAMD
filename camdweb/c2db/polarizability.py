"""Optical polarizability."""

from __future__ import annotations

import typing
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from asr.core import ASRResult, command, option, prepare_result
from camdweb.html import image, table
from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData, SkipPanel
from click import Choice

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
   {images[0]}
   {images[2]}
  </div>
  <div class="col-6">
   {images[1]}
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


class Polarizability(Panel):
    info = INFO

    def get_data(self, material):
        file = material.folder / 'polarization.json'
        if not file.is_file():
            raise SkipPanel

        png = material.folder / 'rpa-pol-x.png'
        if not png.is_file():
            data = json.loads(png.read_text())
            polarizability(data['frequencies'],
                           data['alphax_w'],
                           data['alphay_w'],
                           data['alphaz_w'],
                           material.plasmafrequency_x,
                           material.plasmafrequency_y):

        tab = table(['Properties', ' '],
                    [])

        return PanelData(
            HTML.format(**{f'i{v}': image('rpa-pol-{v}.png')
                           for v in 'xyz'},
                        table=tab),
            title='Optical polarizability')


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


def polarizability(frequencies,
                   alphax_w,
                   alphay_w,
                   alphaz_w,
                   plasmafrequency_x,
                   plasmafrequency_y):
    i2 = abs(frequencies - 50.0).argmin()
    frequencies = frequencies[:i2]
    alphax_w = alphax_w[:i2]
    alphay_w = alphay_w[:i2]
    alphaz_w = alphaz_w[:i2]

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

    fx, fy, fz = (f'rpa-pol-{v}.png' for v in 'xyz')

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
