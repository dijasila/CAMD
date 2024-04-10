from __future__ import annotations
import json
import sys
from typing import Iterable, Generator
from pathlib import Path
import warnings

from camdweb.html import table
from camdweb.html import image
from camdweb.material import Material
from camdweb.c2db.asr_panel import read_result_file as rrf
import numpy as np
from camdweb.panels.panel import Panel

# reference = """\
# S. Kaappa et al. Point group symmetry analysis of the electronic structure
# of bare and protected nanocrystals, J. Phys. Chem. A, 122, 43, 8576 (2018)"""


# panel_description = make_panel_description(
#     """
# Analysis of defect states localized inside the pristine bandgap (energetics and
#  symmetry).
# """,
#     articles=[
#         href(reference, 'https://doi.org/10.1021/acs.jpca.8b07923'),
#     ],
# )

HTML = """
<div class="row">
  <div class="col-6">
    <div id='defect_symmetry' class='defect_symmetry'></div>
    {tbl0}
    {ks_gap_png}
  </div>
  <div class="col-6">
    {tbl1}
    {tbl2}
  </div>
</div>
"""
### Panel class
class DefectSymmetryPanel(Panel):
    title = f'One-electron defect states'
    def get_html(self,
                 material: Material) -> Generator[str, None]:
        
        root = material.folder
        defsym_file = root / 'results-asr.defect_symmetry.json'
        if not defsym_file.is_file():
            return
        
        ks_gap_png = root / 'ks_gap.png'
        tbl0, tbl1, tbl2 = self.plot_and_tables(root, ks_gap_png)
        yield HTML.format(tbl0=tbl0, tbl1=tbl1, tbl2=tbl2,
                          ks_gap_png=image(ks_gap_png, 'One-electon defect states'))
    
    def plot_and_tables(self, root: Path, ks_gap_png: Path):
        data = rrf(root / 'results-asr.defect_symmetry.json')
        e_fermi = rrf(root / 'results-asr.gs.json')['efermi']
        eref = rrf(root / 'results-asr.get_wfs.json')['eref']


        vbm = data['pristine']['kwargs']['data']['vbm']
        cbm = data['pristine']['kwargs']['data']['cbm']
        symmetries_data = [symmetries['kwargs']['data'] 
                           for symmetries in data['symmetries']]

        if symmetries_data[0]['best'] is None:
            warnings.warn("no symmetry analysis present for this defect. "
                      "Only plot gapstates!", UserWarning)
            style = 'state'
        else:
            style = 'symmetry'
            
        state_tables, transition_table = get_symmetry_tables(
            symmetries_data, vbm, cbm, e_fermi, eref, style=style)
        
        plot_gapstates(symmetries_data, cbm, vbm, eref, e_fermi, ks_gap_png)

        return state_tables[0], state_tables[1], transition_table



### Tables and plots for the defect symmetry panel

def get_symmetry_tables(state_results, vbm, cbm, e_fermi, eref, style):
    state_tables = []
    ef = e_fermi - eref # shift to vacuum level

    E_hls = []
    for spin in range(2):
        state_array, rowlabels = get_matrixtable_array(
            state_results, vbm, cbm, spin, style)
        if style == 'symmetry':
            delete = [2]
            columnlabels = ['Symmetry',
                            # 'Spin',
                            'Localization ratio',
                            'Energy']
        elif style == 'state':
            delete = [0, 2, 3]
            columnlabels = [  # 'Spin',
                'Energy']

        N_homo = 0
        N_lumo = 0
        for i in range(len(state_array)):
            if float(state_array[i, 4]) > ef:
                N_lumo += 1

        E_homo = vbm
        E_lumo = cbm
        for i in range(len(state_array)):
            if float(state_array[i, 4]) > ef:
                rowlabels[i] = f'LUMO + {N_lumo - 1}'
                N_lumo = N_lumo - 1
                if N_lumo == 0:
                    rowlabels[i] = 'LUMO'
                    E_lumo = float(state_array[i, 4])
            elif float(state_array[i, 4]) <= ef:
                rowlabels[i] = f'HOMO — {N_homo}'
                if N_homo == 0:
                    rowlabels[i] = 'HOMO'
                    E_homo = float(state_array[i, 4])
                N_homo = N_homo + 1
        E_hl = E_lumo - E_homo
        E_hls.append(E_hl)

        state_array = np.delete(state_array, delete, 1)
        rows = []

        for i in range(len(state_array)):
            if style == 'symmetry':
                rows.append((rowlabels[i],
                             # state_array[i, 0],
                             state_array[i, 1],
                             state_array[i, 2],
                            f'{state_array[i, 3]} eV'))
            elif style == 'state':
                rows.append((rowlabels[i],
                             # state_array[i, 0],
                             f'{state_array[i, 1]} eV'))
        
        state_table = table(
            header=[f'Orbitals in spin channel {spin}',
                    *columnlabels],
            rows=rows)

        #state_table['rows'] = rows
        state_tables.append(state_table)

    transition_table = get_transition_table(E_hls)

    return state_tables, transition_table

def get_transition_table(E_hls):
    """Create table for HOMO-LUMO transition in both spin channels."""

    rows = [[f'Spin {i}', f'{element:.2f} eV']
            for i, element in enumerate(E_hls)]

    transition_table = table(
        header=['Kohn—Sham HOMO—LUMO gap', 'Value'],
        rows=rows)

    return transition_table

def get_number_of_rows(res, spin, vbm, cbm):
    counter = 0
    for i in range(len(res)):
        if (int(res[i]['spin']) == spin
           and res[i]['energy'] < cbm
           and res[i]['energy'] > vbm):
            counter += 1

    return counter

def get_matrixtable_array(state_results, vbm, cbm,
                          spin, style):
    Nrows = get_number_of_rows(state_results, spin, vbm, cbm)
    state_array = np.empty((Nrows, 5), dtype='object')
    rowlabels = []
    spins = []
    energies = []
    symlabels = []
    accuracies = []
    loc_ratios = []
    for i, row in enumerate(state_results):
        rowname = f"{int(state_results[i]['state']):.0f}"
        label = str(state_results[i]['best'])
        labelstr = label.lower()
        splitstr = list(labelstr)
        if len(splitstr) == 2:
            labelstr = f'{splitstr[0]}<sub>{splitstr[1]}</sub>'
        if state_results[i]['energy'] < cbm and state_results[i]['energy'] > vbm:
            if int(state_results[i]['spin']) == spin:
                rowlabels.append(rowname)
                spins.append(f"{int(state_results[i]['spin']):.0f}")
                energies.append(f"{state_results[i]['energy']:.2f}")
                if style == 'symmetry':
                    symlabels.append(labelstr)
                    accuracies.append(f"{state_results[i]['error']:.2f}")
                    loc_ratios.append(f"{state_results[i]['loc_ratio']:.2f}")
    state_array = np.empty((Nrows, 5), dtype='object')
    rowlabels.sort(reverse=True)

    for i in range(Nrows):
        state_array[i, 1] = spins[i]
        if style == 'symmetry':
            state_array[i, 0] = symlabels[i]
            state_array[i, 2] = accuracies[i]
            state_array[i, 3] = loc_ratios[i]
        state_array[i, 4] = energies[i]
    state_array = state_array[state_array[:, -1].argsort()]

    return state_array, rowlabels

## Plotting
class Level:
    """Class to draw a single defect state level in the gap."""

    def __init__(self, energy, spin, deg, off, size=0.05, ax=None):
        self.size = size
        self.energy = energy
        self.ax = ax
        self.spin = spin
        self.deg = deg
        assert deg in [1, 2], ('only degeneracies up to two are '
                               'implemented!')
        self.off = off
        self.relpos = self.get_relative_position(self.spin, self.deg, self.off)

    def get_relative_position(self, spin, deg, off):
        """Set relative position of the level based on spin, degeneracy and offset."""
        xpos_deg = [[2 / 12, 4 / 12], [8 / 12, 10 / 12]]
        xpos_nor = [1 / 4, 3 / 4]
        if deg == 2:
            relpos = xpos_deg[spin][off]
        elif deg == 1:
            relpos = xpos_nor[spin]

        return relpos

    def draw(self):
        """Draw the defect state according to spin and degeneracy."""
        pos = [self.relpos - self.size, self.relpos + self.size]
        self.ax.plot(pos, [self.energy] * 2, '-k')

    def add_occupation(self, length):
        """Draw an arrow if the defect state if occupied."""
        updown = [1, -1][self.spin]
        self.ax.arrow(self.relpos,
                      self.energy - updown * length / 2,
                      0,
                      updown * length,
                      head_width=0.01,
                      head_length=length / 5, fc='C3', ec='C3')

    def add_label(self, label, static=None):
        """Add symmetry label of the irrep of the point group."""
        shift = self.size / 5
        labelcolor = 'C3'
        if static is None:
            labelstr = label.lower()
            splitstr = list(labelstr)
            if len(splitstr) == 2:
                labelstr = f'{splitstr[0]}$_{splitstr[1]}$'
        else:
            labelstr = 'a'

        if (self.off == 0 and self.spin == 0):
            xpos = self.relpos - self.size - shift
            ha = 'right'
        if (self.off == 0 and self.spin == 1):
            xpos = self.relpos + self.size + shift
            ha = 'left'
        if (self.off == 1 and self.spin == 0):
            xpos = self.relpos - self.size - shift
            ha = 'right'
        if (self.off == 1 and self.spin == 1):
            xpos = self.relpos + self.size + shift
            ha = 'left'
        self.ax.text(xpos,
                     self.energy,
                     labelstr,
                     va='center', ha=ha,
                     size=12,
                     color=labelcolor)

def plot_gapstates(symmetries_data, cbm, vbm, eref, e_fermi, fname):
    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()

    # extract pristine data
    gap = cbm - vbm
    ef = e_fermi - eref

    # Draw band edges
    draw_band_edge(vbm, 'vbm', offset=gap / 5, ax=ax)
    draw_band_edge(cbm, 'cbm', offset=gap / 5, ax=ax)

    levelflag = symmetries_data[0]['best'] is not None
    # draw the levels with occupations, and labels for both spins
    for spin in [0, 1]:
        spin_data = get_spin_data(symmetries_data, spin)
        draw_levels_occupations_labels(ax, spin, spin_data, cbm, vbm,
                                       ef, gap, levelflag)

    ax1 = ax.twinx()
    ax.set_xlim(0, 1)
    ax.set_ylim(vbm - gap / 5, cbm + gap / 5)
    ax1.set_ylim(vbm - gap / 5, cbm + gap / 5)
    ax1.plot([0, 1], [ef] * 2, '--k')
    ax1.set_yticks([ef])
    ax1.set_yticklabels([r'$E_\mathrm{F}$'])
    ax.set_xticks([])
    ax.set_ylabel(r'$E-E_\mathrm{vac}$ [eV]')

    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

def get_spin_data(symmetries_data, spin):
    """Create symmetry result only containing entries for one spin channel."""
    spin_data = []
    for sym in symmetries_data:
        if int(sym['spin']) == spin:
            spin_data.append(sym)

    return spin_data

def draw_band_edge(energy, edge, *, offset=2, ax):
    if edge == 'vbm':
        eoffset = energy - offset
        elabel = energy - offset / 2
    elif edge == 'cbm':
        eoffset = energy + offset
        elabel = energy + offset / 2

    ax.plot([0, 1], [energy] * 2, color='black', zorder=1)
    ax.fill_between([0, 1], [energy] * 2, [eoffset] * 2, color='grey', alpha=0.5)
    ax.text(0.5, elabel, edge.upper(), color='w', weight='bold', ha='center',
            va='center', fontsize=12)

def draw_levels_occupations_labels(ax, spin, spin_data, ecbm, evbm, ef,
                                   gap, levelflag):
    """Loop over all states in the gap and plot the levels.

    This function loops over all states in the gap of a given spin
    channel, and draws the states with labels. If there are
    degenerate states, it makes use of the degeneracy_counter, i.e. if two
    degenerate states follow after each other, one of them will be drawn
    on the left side (degoffset=0, degeneracy_counter=0), the degeneracy
    counter will be increased by one and the next degenerate state will be
    drawn on the right side (degoffset=1, degeneracy_counter=1). Since we
    only deal with doubly degenerate states here, the degeneracy counter
    will be set to zero again after drawing the second degenerate state.

    For non degenerate states, i.e. deg = 1, all states will be drawn
    in the middle and the counter logic is not needed.
    """
    # initialize degeneracy counter and offset
    degeneracy_counter = 0
    degoffset = 0
    for sym in spin_data:
        energy = sym['energy']
        is_inside_gap = evbm < energy < ecbm
        if is_inside_gap:
            spin = int(sym['spin'])
            irrep = sym['best']
            # only do drawing left and right if levelflag, i.e.
            # if there is a symmetry analysis to evaluate degeneracies
            if levelflag:
                deg = [1, 2]['E' in irrep]
            else:
                deg = 1
                degoffset = 1
            # draw draw state on the left hand side
            if deg == 2 and degeneracy_counter == 0:
                degoffset = 0
                degeneracy_counter = 1
            # draw state on the right hand side, set counter to zero again
            elif deg == 2 and degeneracy_counter == 1:
                degoffset = 1
                degeneracy_counter = 0
            # intitialize and draw the energy level
            lev = Level(energy, ax=ax, spin=spin, deg=deg,
                        off=degoffset)
            lev.draw()
            # add occupation arrow if level is below E_F
            if energy <= ef:
                lev.add_occupation(length=gap / 15.)
            # draw label based on irrep
            if levelflag:
                static = None
            else:
                static = 'A'
            lev.add_label(irrep, static=static)