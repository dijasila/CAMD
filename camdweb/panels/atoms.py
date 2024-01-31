"""Panel for showing the atomic structure and additional tables.

Content::

  +-----------+-------------------------------+
  | column 1  | repeat-unit-cell button       |
  | ...       | download button (TODO)        |
  |           +-------------------------------+
  |           | ball and stick plotly plot    |
  |           | ...                           |
  |           |                               |
  |           +-------------------------------+
  |           | unit cell: vectors            |
  |           +-------------------------------+
  |           | unit cell: lengths and angles |
  +-----------+-------------------------------+
"""
from __future__ import annotations

import json
import sys
from functools import cache
from math import nan
from typing import Generator

import numpy as np
import plotly
import plotly.graph_objects as go
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.io import read
from ase.neighborlist import neighbor_list
from scipy.spatial import ConvexHull

from camdweb.html import table
from camdweb.material import Material, Materials
from camdweb.panels.panel import Panel

HTML = """
<h4>{formula}</h4>
<div class="row">
  <div class="col-6">
    {column1}
  </div>
  <div class="col-6">
    {column2}
  </div>
</div>
"""


def default_repeat(material):
    cell_cv = material.atoms.cell
    pbc_c = material.atoms.get_pbc()
    if not np.any(pbc_c):
        return 1
    V = np.abs(np.linalg.det(cell_cv[pbc_c][:, pbc_c]))
    return min(4, int(np.round(15 / V**(1 / np.sum(pbc_c)))))


def repeat_options(selected, maximum):
    return "".join([f'<option value="{value}"'
                    f'{" selected" if value == selected else ""}>'
                    f'{value}</option>'
                    for value in range(1, maximum + 1)])


COLUMN2 = """
    <label>Repeat:</label>
    <select onchange="cb(this.value, 'atoms', '{uid}')">
    {options}
    </select>

    &emsp;

    <label>Download:</label>
    <a download="atoms.xyz"
       class="btn btn-info btn-sm" href={uid}/download/xyz>XYZ</a>
    <a download="atoms.cif"
       class="btn btn-info btn-sm" href={uid}/download/cif>CIF</a>
    <a download="atoms.json"
       class="btn btn-info btn-sm" href={uid}/download/json>JSON</a>

    <div id='atoms' class='atoms'></div>
    {axes}

<script type='text/javascript'>
var graphs = {atoms_json};
Plotly.newPlot('atoms', graphs, {{}});
</script>
"""

DIMS = ['', 'length', 'area', 'volume']


class AtomsPanel(Panel):
    title = 'Atoms'

    def __init__(self) -> None:
        self.callbacks = {'atoms': self.plot}
        self.column_names = {}
        self.columns = ['stoichiometry', 'uid']

    def update_data(self, material):
        pbc = material.atoms.pbc
        dims = pbc.sum()
        if dims > 0:
            vol = abs(np.linalg.det(material.atoms.cell[pbc][:, pbc]))
            name = DIMS[dims]
            if name not in self.column_names:
                if dims == 1:
                    unit = 'Å'
                else:
                    unit = f'Å<sup>{dims}</sup>'
                self.column_names[name] = f'{name.title()} [{unit}]'
                self.columns.append(name)
            material.add_column(name, vol)

    def get_html(self,
                 material: Material,
                 materials: Materials) -> Generator[str, None, None]:
        col1 = self.create_column_one(material, materials)
        col2 = self.create_column_two(material, materials)
        yield HTML.format(column1=col1,
                          column2=col2,
                          formula=material['formula'])

    def create_column_one(self,
                          material: Material,
                          materials: Materials) -> str:
        return table(None, materials.table(material, self.columns))

    def create_column_two(self,
                          material: Material,
                          materials: Materials) -> str:
        defrep = default_repeat(material)
        return COLUMN2.format(
            axes=self.axes(material),
            options=repeat_options(defrep, 4),
            uid=material.uid,
            atoms_json=self.plot(material, defrep))

    def axes(self, material: Material) -> str:
        atoms = material.atoms
        tbl1 = table(
            ['Axis', 'x [Å]', 'y [Å]', 'z [Å]', 'Periodic'],
            [[i + 1, *[f'{x:.3f}' for x in axis], 'Yes' if p else 'No']
             for i, (axis, p) in enumerate(zip(atoms.cell, atoms.pbc))])
        C = atoms.cell.cellpar()
        tbl2 = table(
            None,
            [['Lengths [Å]', *[f'{x:.3f}' for x in C[:3]]],
             ['Angles [°]', *[f'{x:.3f}' for x in C[3:]]]])
        return tbl1 + tbl2

    def plot(self, material: Material, repeat: int = 1) -> str:
        assert repeat < 5, 'DOS!'
        atoms = material.atoms
        unitcell = atoms.cell.copy()
        atoms2 = atoms.copy()
        show_magmoms = True
        try:
            atoms2.set_initial_magnetic_moments(atoms.calc.results['magmoms'])
        except (AttributeError, KeyError):
            show_magmoms = False
        atoms2 = atoms2.repeat([repeat if p else 1 for p in atoms.pbc])
        fig = plot_atoms(atoms2, unitcell, show_magmoms=show_magmoms)
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


# The 12 edges of a cube:
UNITCELL = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1],
            [nan, nan, nan], [1, 0, 0], [1, 0, 1],
            [nan, nan, nan], [0, 1, 0], [0, 1, 1],
            [nan, nan, nan], [1, 1, 0], [1, 1, 1]]


def get_bonds(atoms):
    i, j, D, S, d = neighbor_list('ijDSd', atoms,
                                  cutoff=covalent_radii[atoms.numbers] * 1.2)

    # For bonds inside the cell (i.e., cell displacement S==(0, 0, 0)),
    # bonds are double-counted.  This mask un-doublecounts them:
    doublecount_mask = (i < j) | S.any(1)
    i = i[doublecount_mask]
    j = j[doublecount_mask]
    D = D[doublecount_mask]
    S = S[doublecount_mask]
    d = d[doublecount_mask]
    return i, j, D, S, d


def plot_atoms(atoms: Atoms,
               unitcell: np.ndarray | None = None,
               show_magmoms: bool = False) -> go.Figure:
    """Ball and stick plotly figure."""
    data = []
    # Atoms:
    points = SPHERE_POINTS
    triangles = triangulate_sphere()
    i, j, k = triangles.T
    for Z, symbol, xyz, magmom in zip(atoms.numbers,
                                      atoms.symbols,
                                      atoms.positions,
                                      atoms.get_initial_magnetic_moments()):
        x, y, z = (points * covalent_radii[Z] * 0.5 + xyz).T
        magmoms = f'<i>Magmom</i> {magmom:.2f}' if show_magmoms else ''
        hovertemplate = (
            f'<b>{symbol}</b><br><i>Coord.</i>'
            f'{xyz[0]:.2f} {xyz[1]:.2f} {xyz[2]:.2f}<br>{magmoms}')

        mesh = go.Mesh3d(x=x, y=y, z=z,
                         i=i, j=j, k=k, opacity=1,
                         color=color(Z),
                         text=symbol,
                         name='',
                         hovertemplate=hovertemplate)
        data.append(mesh)

    # Bonds:
    i, j, D, S, d = get_bonds(atoms)
    p = atoms.positions

    # Only add hover text to bond center
    text = []
    for L in d:
        bond_txt = f'{L:.2f} Å'
        text += ['', bond_txt, '', '']

    xyz = np.empty((3, len(i) * 4))

    # bond consists of 3 actual points (0, 1 and 2) and final nan value (3)
    # to break the line, and this is repeated in sequence for each bond
    xyz[:, 0::4] = p[i].T
    xyz[:, 1::4] = (p[i] + D * 0.5).T
    xyz[:, 2::4] = (p[i] + D).T
    # If bond exits supercell, shorten it by half
    xyz[:, 2::4][:, S.any(1)] = np.nan
    xyz[:, 3::4] = np.nan

    x, y, z = xyz
    data.append(go.Scatter3d(x=x, y=y, z=z, text=text, mode='lines',
                             hovertemplate='%{text}',
                             name='',
                             line=dict(color='grey', width=5),
                             showlegend=False))

    # Unit cell:
    if unitcell is None:
        unitcell = atoms.cell
    x, y, z = (UNITCELL @ unitcell).T
    data.append(go.Scatter3d(x=x, y=y, z=z, mode='lines',
                             hovertemplate='',
                             hoverinfo='skip',
                             name='',
                             line=dict(color='#fa9fb5', width=10),
                             showlegend=False))

    fig = go.Figure(data=data)
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.update_xaxes(showgrid=False)
    fig.update_layout(template='simple_white')
    fig.update_scenes(aspectmode='data',
                      camera=dict(projection=dict(type='orthographic')))
    return fig


# 50 Lebedev quadrature points:
SPHERE_POINTS = np.array(
    [[-1.0, 0.0, 0.0],
     [1.0, 0.0, 0.0],
     [0.0, -1.0, 0.0],
     [0.0, 1.0, 0.0],
     [0.0, 0.0, -1.0],
     [0.0, 0.0, 1.0],
     [0.0, -0.70710678118654757, -0.70710678118654757],
     [0.0, -0.70710678118654757, 0.70710678118654757],
     [0.0, 0.70710678118654757, -0.70710678118654757],
     [0.0, 0.70710678118654757, 0.70710678118654757],
     [-0.70710678118654757, 0.0, -0.70710678118654757],
     [0.70710678118654757, 0.0, -0.70710678118654757],
     [-0.70710678118654757, 0.0, 0.70710678118654757],
     [0.70710678118654757, 0.0, 0.70710678118654757],
     [-0.70710678118654757, -0.70710678118654757, 0.0],
     [-0.70710678118654757, 0.70710678118654757, 0.0],
     [0.70710678118654757, -0.70710678118654757, 0.0],
     [0.70710678118654757, 0.70710678118654757, 0.0],
     [-0.57735026918962573, -0.57735026918962573, -0.57735026918962573],
     [-0.57735026918962573, -0.57735026918962573, 0.57735026918962573],
     [-0.57735026918962573, 0.57735026918962573, -0.57735026918962573],
     [-0.57735026918962573, 0.57735026918962573, 0.57735026918962573],
     [0.57735026918962573, -0.57735026918962573, -0.57735026918962573],
     [0.57735026918962573, -0.57735026918962573, 0.57735026918962573],
     [0.57735026918962573, 0.57735026918962573, -0.57735026918962573],
     [0.57735026918962573, 0.57735026918962573, 0.57735026918962573],
     [-0.90453403373329089, -0.30151134457776357, -0.30151134457776357],
     [-0.90453403373329089, -0.30151134457776357, 0.30151134457776357],
     [-0.90453403373329089, 0.30151134457776357, -0.30151134457776357],
     [-0.90453403373329089, 0.30151134457776357, 0.30151134457776357],
     [0.90453403373329089, -0.30151134457776357, -0.30151134457776357],
     [0.90453403373329089, -0.30151134457776357, 0.30151134457776357],
     [0.90453403373329089, 0.30151134457776357, -0.30151134457776357],
     [0.90453403373329089, 0.30151134457776357, 0.30151134457776357],
     [-0.30151134457776357, -0.90453403373329089, -0.30151134457776357],
     [0.30151134457776357, -0.90453403373329089, -0.30151134457776357],
     [-0.30151134457776357, -0.90453403373329089, 0.30151134457776357],
     [0.30151134457776357, -0.90453403373329089, 0.30151134457776357],
     [-0.30151134457776357, 0.90453403373329089, -0.30151134457776357],
     [0.30151134457776357, 0.90453403373329089, -0.30151134457776357],
     [-0.30151134457776357, 0.90453403373329089, 0.30151134457776357],
     [0.30151134457776357, 0.90453403373329089, 0.30151134457776357],
     [-0.30151134457776357, -0.30151134457776357, -0.90453403373329089],
     [-0.30151134457776357, 0.30151134457776357, -0.90453403373329089],
     [0.30151134457776357, -0.30151134457776357, -0.90453403373329089],
     [0.30151134457776357, 0.30151134457776357, -0.90453403373329089],
     [-0.30151134457776357, -0.30151134457776357, 0.90453403373329089],
     [-0.30151134457776357, 0.30151134457776357, 0.90453403373329089],
     [0.30151134457776357, -0.30151134457776357, 0.90453403373329089],
     [0.30151134457776357, 0.30151134457776357, 0.90453403373329089]])


@cache
def triangulate_sphere() -> np.ndarray:
    hull = ConvexHull(SPHERE_POINTS)
    tri_tv = hull.simplices.copy()

    # Make sure surface normals are pointing out from the sphere surface
    tri_tvc = SPHERE_POINTS[tri_tv]
    n_tc = np.cross(tri_tvc[:, 1, :] - tri_tvc[:, 0, :],
                    tri_tvc[:, 2, :] - tri_tvc[:, 0, :])
    flip_n = np.sum(n_tc * tri_tvc[:, 0, :], axis=1) < 0
    tri_tv[flip_n, 1], tri_tv[flip_n, 2] = tri_tv[flip_n, 2], tri_tv[flip_n, 1]

    return tri_tv


@cache
def color(Z: int) -> str:
    a, b, c = (jmol_colors[Z] * 255).astype(int)
    return f'rgb({a},{b},{c})'


if __name__ == '__main__':
    atoms = read(sys.argv[1])
    assert isinstance(atoms, Atoms)
    plot_atoms(atoms).show()
