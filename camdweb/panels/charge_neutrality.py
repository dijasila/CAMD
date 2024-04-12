from __future__ import annotations
import json
import sys
import plotly
from typing import Iterable

from camdweb.html import table
from camdweb.material import Material
from camdweb.c2db.asr_panel import read_result_file
from camdweb.panels.panel import Panel, PanelData
import numpy as np
import plotly.graph_objects as go

# panel_description = make_panel_description(
#     """
# Equilibrium defect energetics evaluated by solving E<sub>F</sub> self-consistently
# until charge neutrality is achieved.
# """,
#     articles=[
#         href("""J. Buckeridge, Equilibrium point defect and charge carrier
#  concentrations in a meterial determined through calculation of the self-consistent
#  Fermi energy, Comp. Phys. Comm. 244 329 (2019)""",
#              'https://doi.org/10.1016/j.cpc.2019.06.017'),
#     ],
# )

HTML = """
<div class="row">
  <div class="col-6">
    <div id='{div_id}' class='charge_neutrality'></div>
    {tbl0}
  </div>
  <div class="col-6">
    {tbl1}
  </div>
</div>
"""

SCRIPT = """
<script type='text/javascript'>
{graphs}
</script>
"""

class ChargeNeutralitySuperpanel(Panel):
    def __init__(self) -> None:
        super().__init__()
        #self.subpanels: list[Panel] = []

    def add_subpanels(self, material: Material):
        root = material.folder.parent.parent
        cn_file = root / 'pristine_sc' / 'results-asr.charge_neutrality.json'
        result = read_result_file(cn_file)

        subpanels = list()

        # Iterate over each "condition" in "scresults"
        for i, scresult in enumerate(result['scresults']):
            scresult = scresult['kwargs']['data']
            condition = scresult['condition']

            subpanels.append(
                ChargeNeutralityPanel(result, scresult, condition).get_data(material))
        return subpanels

    def get_data(self, material: Material) -> PanelData:
        subpanels = self.add_subpanels(material)
        return PanelData(html='', title='Equilibrium energetics: All defects', subpanels=subpanels)


class ChargeNeutralityPanel(Panel):
    def __init__(self, result, scresult, condition) -> None:
        super().__init__()
        self.result = result
        self.scresult = scresult
        self.title = condition

    def get_data(self,
                 material: Material) -> PanelData:

        html=''
        graphs_js = ''

        result = self.result
        scresult = self.scresult
        title = self.title

        charge_neutrality_fig, tbl0, tbl1 = make_figure_and_tables(scresult,
                                                                    result,
                                                                    verbose=False)

        charge_neutrality_json = json.dumps(charge_neutrality_fig,
                                            cls=plotly.utils.PlotlyJSONEncoder)
        html += f'<h3>{title}</h3>'
        html += HTML.format(div_id=f'charge_neutrality_{title}', tbl0=tbl0, tbl1=tbl1)

        graphs_js = f"var graph{title} = {charge_neutrality_json};\n"
        graphs_js += f"Plotly.newPlot('charge_neutrality_{title}', graph{title}, {{}});\n"

        html += SCRIPT.format(graphs=graphs_js) + '<br>'

        return PanelData(html,
                    title=title,
                    script=graphs_js)

# Make equilibrium defects energetics plot (from asr.charge_neutrality), //
# an overview table for SCF results, and a table for each defect's equilibrium concentration
def make_figure_and_tables(scresult: dict[str, tuple[dict[str, int],
                                                 float,
                                                 str]],
                            result: dict[str, tuple[dict[str, int],
                                                 float,
                                                 str]],
                           verbose: bool = True) -> tuple[str, str, str]:

    unit = result['conc_unit']
    assert isinstance(unit, str)
    unitstring = f"cm<sup>{unit.split('^')[-1]}</sup>"

    """Make charge neutrality figure and tables."""
    condition = scresult['condition']
    plotname = f'neutrality-{condition}.png'
    fig = plot_formation_scf(scresult, result, plotname)

    tbl0 = get_overview_table(scresult, result, unitstring)

    tbl1 = [] # this is list of tables [[(header1, rows1)],[(header2, rows2)], ...]

    conc_rows = []
    conc_headers = []
    for defect_concentration in scresult['defect_concentrations']:
        element = defect_concentration['kwargs']['data']
        if isinstance(element, str):
            element = json.loads(element)
        conc_row, def_name, def_type = get_conc_table(result, element, unitstring)
        conc_header = ['Charge state', f'Eq. concentrations of '
                       f'{def_type}<sub>{def_name}</sub> [{unitstring}]']
        conc_rows.append(conc_row)
        conc_headers.append(conc_header)

    scf_temp = result['temperature'] # scf temperature

    html0 = table([f'Equilibrium properties @ {int(scf_temp):d} K', 'Value'],
                  [(item1, item2) for item1, item2 in tbl0])

    # Create many tables for each conc_header \\
    # And iterate over conc_headers, with rows conc_rows[rows] \\
    # to create tables like table([header, rows])
    html1 = ''
    for header, rows in zip(conc_headers, conc_rows):
        html1 += table(header, rows)

    return fig, html0, html1

def get_overview_table(scresult, result, unitstring):
    ef = scresult['efermi_sc']
    gap = result['gap']
    if ef < (gap / 4.):
        dopability = '<b style="color:red;">p-type</b>'
    elif ef > (3 * gap / 4.):
        dopability = '<b style="color:blue;">n-type</b>'
    else:
        dopability = 'intrinsic'

    # get strength of p-/n-type dopability
    if ef < 0:
        ptype_val = '100+'
        ntype_val = '0'
    elif ef > gap:
        ptype_val = '0'
        ntype_val = '100+'
    else:
        ptype_val = int((1 - ef / gap) * 100)
        ntype_val = int((100 - ptype_val))

    pn_strength = f'{ptype_val:3}% / {ntype_val:3}%'
    pn = 'Strength of p-/n-type dopability in percent'

    is_dopable = 'Intrinsic doping type'

    scf_fermi = 'Fermi level position'

    scf_overview = []
    scf_overview.append((is_dopable, dopability))
    scf_overview.append((scf_fermi, f'{ef:.2f} eV'))
    scf_overview.append((pn, pn_strength))
    if scresult['n0'] > 1e-5:
        n0 = scresult['n0']
    else:
        n0 = 0
    scf_overview.append(
        ('Electron carrier concentration',
          f'{n0:.1e} {unitstring}'))
    if scresult['p0'] > 1e-5:
        p0 = scresult['p0']
    else:
        p0 = 0
    scf_overview.append(
        ('Hole carrier concentration',
          f'{p0:.1e} {unitstring}'))

    return scf_overview

def get_conc_table(result, element, unitstring):
    name = element['defect_name']
    def_type = name.split('_')[0]
    if def_type == 'v':
        def_type = 'V'
    def_name = name.split('_')[1]
    scf_conc_table = []
    for altel in element['concentrations']:
        if altel[0] > 1e1:
            scf_conc_table.append(
                (f'<b>Charge {altel[1]:1d}</b>',
                  f'<b>{altel[0]:.1e}</b>'))
        else:
            scf_conc_table.append(
                (f'Charge {altel[1]:1d}',
                  f'{altel[0]:.1e}'))

    return scf_conc_table, def_name, def_type

def plot_formation_scf(scresult, result, fname) -> go.Figure:
    """Plot formation energy diagram and SC Fermi level wrt. VBM."""
    gap = result['gap']
    comparison = fname.split('neutrality-')[-1].split('.png')[0]
    fig = go.Figure()

    if comparison == scresult['condition']:
        ef = scresult['efermi_sc']
        for i, defect in enumerate(scresult['defect_concentrations']):
            defect = defect['kwargs']['data']
            name = defect['defect_name']
            def_type = name.split('_')[0]
            def_name = name.split('_')[-1]
            if def_type == 'v':
                def_type = 'V'
            namestring = f"{def_type}<sub>{def_name}</sub>"
            array = np.zeros((len(defect['concentrations']), 2))
            for num, conc_tuple in enumerate(defect['concentrations']):
                q = conc_tuple[1]
                eform = conc_tuple[2]
                array[num, 0] = eform + q * (-ef)
                array[num, 1] = q
            array = array[array[:, 1].argsort()[::-1]]
            plot_lowest_lying(fig, array, gap, name=namestring)
        draw_band_edges(fig, gap)
        set_limits(fig, gap)
        draw_ef(fig, ef)
        set_labels_and_legend(fig, comparison)
    return fig

# Helper functions for plotting
def set_labels_and_legend(fig, title):
    fig.update_layout(
        xaxis=dict(title='E<sub>F</sub> - E<sub>VBM</sub> [eV]'),
        yaxis=dict(title='E<sup>f</sup> [eV]'),
        title=title,
        legend=dict(x=0.5, y=1.1, orientation='h',
                    xanchor='center', yanchor='bottom',
                    bgcolor='rgba(0,0,0,0)'))

def draw_ef(fig, ef):
    fig.add_shape(
        type="line", x0=ef, y0=0, x1=ef, y1=1,
        yref="paper", # y in normalized coordinates with respect to the paper
        line=dict(color="red", dash="dot"),
        name="E<sub>F</sub><sup>sc</sup>",
    )
    fig.update_yaxes(autorange=True)

def set_limits(fig, gap):
    fig.update_layout(
        xaxis=dict(range=[0 - gap / 10., gap + gap / 10.])
    )

def get_min_el(array):
    index = np.argmin(array[:, 0])
    return index

def get_crossing_point(y1, y2, q1, q2):
    """
    f1 = y1 + x * q1
    f2 = y2 + x * q2

    x * (q1 - q2) = y2 - y1
    x = (y2 - y1) / (q1 - q2)
    """
    return (y2 - y1) / float(q1 - q2)

def clean_array(array):
    index = get_min_el(array)

    return array[index:, :]

def get_y(x, array, index):
    q = array[index, 1]

    return q * x + array[index, 0]


def get_last_element(array, x_axis, y_axis, gap):
    y_cbms = []
    for i in range(len(array)):
        q = array[i, 1]
        eform = array[i, 0]
        y_cbms.append(q * gap + eform)

    x_axis.append(gap)
    y_axis.append(min(y_cbms))

    return x_axis, y_axis


def get_line_segment(array, index, x_axis, y_axis, gap):
    xs = []
    ys = []
    for i in range(len(array)):
        if i > index:
            y1 = array[index, 0]
            q1 = array[index, 1]
            y2 = array[i, 0]
            q2 = array[i, 1]
            crossing = get_crossing_point(y1, y2, q1, q2)
            xs.append(crossing)
            ys.append(q1 * crossing + y1)
        else:
            crossing = 1000
            xs.append(gap + 10)
            ys.append(crossing)
    min_index = index + 1
    for i, x in enumerate(xs):
        q1 = array[index, 1]
        y1 = array[index, 0]
        if x == min(xs) and x > 0 and x < gap:
            min_index = i
            x_axis.append(xs[min_index])
            y_axis.append(q1 * xs[min_index] + y1)

    return min_index, x_axis, y_axis

def plot_lowest_lying(fig, array_in, gap, name):
    array_tmp = array_in.copy()
    array_tmp = clean_array(array_tmp)
    xs = [0]
    ys = [array_tmp[0, 0]]
    index, xs, ys = get_line_segment(array_tmp, 0, xs, ys, gap)
    for i in range(len(array_tmp)):
        index, xs, ys = get_line_segment(array_tmp, index, xs, ys, gap)
        if index == len(array_tmp):
            break
    xs, ys = get_last_element(array_tmp, xs, ys, gap)

    fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name=name))
    fig.update_layout(xaxis_title="E<sub>F</sub> [eV]")
    fig.write_html('plt_lowest_lying.html')

def draw_band_edges(fig, gap):
    fig = go.Figure()
    fig.add_shape(type="line", x0=0, y0=0, x1=0, y1=1, line=dict(color='black'))
    fig.add_shape(type="line", x0=gap, y0=0, x1=gap, y1=1, line=dict(color='black'))
    fig.add_shape(type="rect", x0=-100, y0=0, x1=0, y1=1, fillcolor='grey', opacity=0.5)
    fig.add_shape(type="rect", x0=gap, y0=0, x1=100, y1=1, fillcolor='grey', opacity=0.5)
