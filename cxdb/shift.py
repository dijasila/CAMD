import numpy as np
from cxdb.panel import Panel, creates
from cxdb.material import Material

HTML = """
% for
<img alt="DOS for {uid}" src="/png/{uid}/dos.png" />
"""


class ShiftPanel(Panel):
    title = 'Shift current spectrum (RPA)'

    def get_html(self,
                 material: Material,
                 column_names: dict[str, str]) -> tuple[str, str]:
        self.make_figures(material)
        return (HTML.format(uid=material.uid), '')

    @creates('shift1.png')
    def make_figures(self, material):
        data = ...  # row.data.get('results-asr.shift.json')

        # Make the table
        sym_chi = data.get('symm')
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
                relation_new = ''  # '\n'.join(wrap(relation, 50))
            table.append((pol, relation_new))
        """
        opt = {'type': 'table',
               'header': ['Element', 'Relations'],
               'rows': table}

        # Make the figure list
        npan = len(sym_chi) - 1
        files = ['shift{}.png'.format(ii + 1) for ii in range(npan)]
        cols = [[fig(f'shift{2 * ii + 1}.png'),
                 fig(f'shift{2 * ii + 2}.png')] for ii in range(int(npan / 2))]
        if npan % 2 == 0:
            cols.append([opt, None])
        else:
            cols.append([fig(f'shift{npan}.png'), opt])
        """


def plot_shift(row, *filename):
    import matplotlib.pyplot as plt
    import os
    from pathlib import Path
    from textwrap import wrap

    # Read the data from the disk
    data = row.data.get('results-asr.shift.json')
    gap = row.get('gap_dir_nosoc')
    atoms = row.toatoms()
    pbc = atoms.pbc.tolist()
    nd = np.sum(pbc)
    if data is None:
        return

    # Remove the files if it is already exist
    for fname in filename:
        if (Path(fname).is_file()):
            os.remove(fname)

    # Plot the data and add the axis labels
    sym_chi = data['symm']
    if len(sym_chi) == 1:
        raise ValueError  # CentroSymmetric
    sigma = data['sigma']

    if not sigma:
        return
    w_l = data['freqs']
    fileind = 0
    axes = []

    for pol in sorted(sigma.keys()):
        # Make the axis and add y=0 axis
        shift_l = sigma[pol]
        ax = plt.figure().add_subplot(111)
        ax.axhline(y=0, color='k')

        # Add the bandgap
        if gap is not None:
            ax.axvline(x=gap, color='k', ls='--')

        # Plot the data
        ax.plot(w_l, np.real(shift_l), '-', c='C0',)

        # Set the axis limit
        ax.set_xlim(0, np.max(w_l))
        relation = sym_chi.get(pol)
        if not (relation is None):
            figtitle = '$' + '$\n$'.join(wrap(relation, 40)) + '$'
            ax.set_title(figtitle)
        ax.set_xlabel(r'Energy [eV]')
        polstr = f'{pol}'
        if nd == 2:
            ax.set_ylabel(r'$\sigma^{(2)}_{' + polstr + r'}$ [nm$\mu$A/V$^2$]')
        else:
            ax.set_ylabel(r'$\sigma^{(2)}_{' + polstr + r'} [$\mu$A/V$^2$]')
        ax.ticklabel_format(axis='both', style='plain', scilimits=(-2, 2))

        # Remove the extra space and save the figure
        plt.tight_layout()
        plt.savefig(filename[fileind])
        fileind += 1
        axes.append(ax)
        plt.close()
