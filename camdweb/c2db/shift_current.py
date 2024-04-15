from pathlib import Path
from textwrap import wrap
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from camdweb.c2db.asr_panel import read_result_file
from camdweb.html import image, table
from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData

HTML = """
<div class="row">
  <div class="col-6">
   {col1}
  </div>
  <div class="col-6">
   {col2}
  </div>
</div>
"""


class ShiftCurrentPanel(Panel):
    datafiles = ['results-asr.shift.json']
    info = """The frequency-dependent shift current is a second order
    optical response tensor, relating a DC current to the square of an
    applied electric field of frequency omega.  The shift current tensor
    is of rank 3 corresponding to the polarization direction of the
    induced current and the two electric field vectors, respectively.  The
    shift current spectra are obtained from second-order perturbation
    theory in the long wave length limit (q=0) without spin-orbit
    interactions.  The calculations include a line-shape broadening of
    50meV, and transitions of near degenerate bands are ignored when the
    energy difference is below 10meV.

    Relevant articles:

    M. O. Sauer, et al.,
    Shift current photovoltaic efficiency in 2D materials.
    npj Comput. Mater. 9, 35 (2023).
    """

    def get_data(self,
                 material: Material) -> PanelData:
        result_file = material.folder / self.datafiles[0]
        data = read_result_file(result_file)
        gap = material.gap_dir_nosoc
        table_rows, figures = make_figures(data, gap, material.folder)
        tbl = table(['Element', 'Relations'], table_rows)
        col1 = [image(path, alt='Shift-current') for path in figures[::2]]
        col2 = [image(path, alt='Shift-current') for path in figures[1::2]]
        if len(figures) % 2 == 0:
            col1.append(tbl)
        else:
            col2.append(tbl)
        return PanelData(HTML.format(col1='\n'.join(col1),
                                     col2='\n'.join(col2)),
                         title='Shift current spectrum (RPA)')


def make_figures(data: dict[str, Any],
                 gap: float,
                 folder: Path) -> tuple[list[tuple[str, str]],
                                        list[Path]]:
    # Make the table
    sym_chi = data['symm']
    table = []
    for pol in sorted(sym_chi.keys()):
        relation = sym_chi[pol]
        if pol == 'zero':
            if relation != '':
                pol = 'Others'
                relation = '0=' + relation
            else:
                continue

        if (len(relation) == 3):
            relation_new = ''
        else:
            # relation_new = '$'+'$\n$'.join(wrap(relation, 40))+'$'
            relation_new = '\n'.join(wrap(relation, 50))
        table.append((pol, relation_new))

    # Make the figure list
    npan = len(sym_chi) - 1
    files = [folder / f'shift{ii + 1}.png' for ii in range(npan)]

    if not files[0].is_file():
        plot_shift(data, gap, files, nd=2)

    return table, files


def plot_shift(data, gap, filenames, nd=2):
    # Plot the data and add the axis labels
    sym_chi = data['symm']
    assert len(sym_chi) != 1, sym_chi  # CentroSymmetric
    sigma = data['sigma']
    if not sigma:
        return  # pragma: no cover
    w_l = data['freqs']

    axes = []

    for filename, pol in zip(filenames, sorted(sigma.keys())):
        # Make the axis and add y=0 axis
        shift_l = sigma[pol]
        ax = plt.figure().add_subplot(111)
        ax.axhline(y=0, color='k')

        # Add the bandgap
        ax.axvline(x=gap, color='k', ls='--')

        # Plot the data
        ax.plot(w_l, np.real(shift_l), '-', c='C0',)

        # Set the axis limit
        ax.set_xlim(0, np.max(w_l))
        relation = sym_chi.get(pol)
        figtitle = '$' + '$\n$'.join(wrap(relation, 40)) + '$'
        ax.set_title(figtitle)
        ax.set_xlabel(r'Energy [eV]')
        polstr = f'{pol}'
        if nd == 2:
            ax.set_ylabel(r'$\sigma^{(2)}_{' + polstr + r'}$ [nm$\mu$A/V$^2$]')
        else:  # pragma: no cover
            ax.set_ylabel(r'$\sigma^{(2)}_{' + polstr + r'} [$\mu$A/V$^2$]')
        ax.ticklabel_format(axis='both', style='plain', scilimits=(-2, 2))

        # Remove the extra space and save the figure
        plt.tight_layout()
        plt.savefig(filename)
        axes.append(ax)
        plt.close()
