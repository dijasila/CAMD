from __future__ import annotations
import json
import sys
from typing import Iterable, Generator
from pathlib import Path

from camdweb.html import table
from camdweb.html import image
from camdweb.material import Material
from camdweb.c2db.asr_panel import read_result_file
import numpy as np
from camdweb.panels.panel import Panel, PanelData, SkipPanel

# panel_description = make_panel_description(
#     """
# Analysis of the thermodynamic stability of the defect using Slater-Janak
#  transition state theory.
# """,
#     articles=[
#         href("""M. Pandey et al. Defect-tolerant monolayer transition metal
# dichalcogenides, Nano Letters, 16 (4) 2234 (2016)""",
#              'https://doi.org/10.1021/acs.nanolett.5b04513'),
#     ],
# )

HTML = """
<div class="row">
  <div class="col-6">
    <div id='slater-janak' class='slater-janak'></div>
    {sj_png}
    {tbl0}
  </div>
  <div class="col-6">
    {tbl1}
  </div>
</div>
"""


class SlaterJanakPanel(Panel):
    def __init__(self) -> None:
        super().__init__()

    def get_data(self,
                 material: Material) -> PanelData:

        root = material.folder.parent
        sj_file = root / 'charge_0' / 'results-asr.sj_analyze.json'
        if not sj_file.is_file():
            raise SkipPanel

        sj_png = root / 'charge_0' / 'sj_transitions.png'
        tbl0, tbl1 = self.plot_and_tables(sj_file, sj_png)

        html = HTML.format(tbl0=tbl0, tbl1=tbl1,
                          sj_png=image(sj_png, 'Slater-Janak charge transitions'))

        return PanelData(html,
                    title = 'Formation energies and Slater-Janak charge transitions',
                    script='')

    def plot_and_tables(self, sj_file: Path, sj_png: Path):
        defect_name = sj_file.parent.parent.name
        data = read_result_file(sj_file)
        vbm = data['pristine']['kwargs']['data']['vbm']

        transitions_data = []
        for transition in data['transitions']:
          transition_data = transition['kwargs']['data']
          transition_name = transition_data['transition_name']
          transition_values = transition_data['transition_values']['kwargs']['data']
          transitions_data.append((transition_name, transition_values))


        ## tbl0 rows, sorted by transition energy
        rows0 = [(f"{defect_name} ({name})",
                 f"{(values['transition'] - values['evac'] - vbm):.2f}",
                 f"{values['erelax']:.2f}")
                for name, values in transitions_data]
        rows0.sort(key=lambda row: float(row[1]))

        tbl0 = table(
            header=['Transition', 'Transition energy', 'Relaxation correction'],
            rows=rows0)

        tbl1 = table(
            header=['Defect formation energy', 'Value'],
            rows=[(f"{defect_name} (q = {row[1]} @ VBM)", f"{row[0]:.2f}")
                  for row in data['eform']])

        # if not sj_png.is_file():
        plot_charge_transitions(data, transitions_data, sj_png)

        return tbl0, tbl1


def plot_charge_transitions(data, transitions_data, sj_png: Path):
    """Plot calculated CTL along with the pristine bandgap."""
    import matplotlib.pyplot as plt

    colors = {'0': 'C0',
              '1': 'C1',
              '2': 'C2',
              '3': 'C3',
              '-1': 'C4',
              '-2': 'C5',
              '-3': 'C6',
              '-4': 'C7',
              '4': 'C8'}

    vbm = data['pristine']['kwargs']['data']['vbm']
    cbm = data['pristine']['kwargs']['data']['cbm']
    gap = abs(cbm - vbm)

    plt.xlim(-1, 1)
    plt.ylim(-0.2 * gap, gap + 0.2 * gap)
    plt.xticks([], [])

    plt.axhspan(-5, 0, color='grey', alpha=0.5)
    plt.axhspan(gap, gap + 5, color='grey', alpha=0.5)
    plt.axhline(0, color='black', linestyle='solid')
    plt.axhline(gap, color='black', linestyle='solid')
    plt.text(0, -0.1 * gap, 'VBM', color='white',
             ha='center', va='center', weight='bold')
    plt.text(0, gap + 0.1 * gap, 'CBM', color='white',
             ha='center', va='center', weight='bold')

    for name, trans_data in transitions_data:
        q = int(name.split('/')[-1])
        q_new = int(name.split('/')[0])
        if q > 0 and q <= 3:
            y = (trans_data['transition']
                 - trans_data['erelax']   # Changed sign /ks
                 - trans_data['evac'])
            color1 = colors[str(q)]
            color2 = colors[str(q_new)]
        elif q < 0 and q >= -3:
            y = (trans_data['transition']
                 + trans_data['erelax']   # Changed sign /ks
                 - trans_data['evac'])
            color1 = colors[str(q)]
            color2 = colors[str(q_new)]

        if -3 <= q <= 3:
            if y <= (cbm + 0.2 * gap) and y >= (vbm - 0.2 * gap):
              plt.plot(np.linspace(-0.9, 0.5, 20), 20 * [y - vbm], label=name,
                       color=color1, mec=color2, mfc=color2, marker='s', markersize=3)
        else:
            continue
    plt.legend(loc='center right')
    plt.ylabel(r'$E - E_{\mathrm{VBM}}$ [eV]')
    plt.yticks()
    plt.tight_layout()
    plt.savefig(sj_png, bbox_inches='tight')
    plt.close()