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

import numpy as np
import plotly
import plotly.graph_objects as go
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.io import read
from ase.neighborlist import neighbor_list
from scipy.spatial import ConvexHull

from cxdb.material import Material, Materials
from cxdb.panels.panel import Panel
from cxdb.utils import table

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

COLUMN2 = """
    <label>Repeat:</label>
    <select onchange="cb(this.value, 'atoms', '{uid}')">
      <option value="1">1</option>
      <option value="2">2</option>
      <option value="3" selected>3</option>
    </select>
    <div id='atoms' class='atoms'></div>
    {axes}
"""

FOOTER = """
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
                 materials: Materials) -> tuple[str, str]:
        col1, foot1 = self.create_column_one(material, materials)
        col2, foot2 = self.create_column_two(material, materials)
        return (HTML.format(column1=col1,
                            column2=col2,
                            formula=material['formula']),
                foot1 + foot2)

    def create_column_one(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        return (table(None, materials.table(material, self.columns)), '')

    def create_column_two(self,
                          material: Material,
                          materials: Materials) -> tuple[str, str]:
        return (
            COLUMN2.format(
                axes=self.axes(material),
                uid=material.uid),
            FOOTER.format(atoms_json=self.plot(material, 3)))

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
        atoms = material.atoms.repeat([repeat if p else 1 for p in atoms.pbc])
        fig = plot_atoms(atoms, unitcell)
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


# The 12 edges of a cube:
UNITCELL = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1],
            [nan, nan, nan], [1, 0, 0], [1, 0, 1],
            [nan, nan, nan], [0, 1, 0], [0, 1, 1],
            [nan, nan, nan], [1, 1, 0], [1, 1, 1]]


def get_bonds(atoms):
    i, j, D, S = neighbor_list('ijDS', atoms,
                               cutoff=covalent_radii[atoms.numbers] * 1.2)

    # For bonds inside the cell (i.e., cell displacement S==(0, 0, 0)),
    # bonds are double-counted.  This mask un-doublecounts them:
    doublecount_mask = (i < j) | S.any(1)
    i = i[doublecount_mask]
    j = j[doublecount_mask]
    D = D[doublecount_mask]
    S = S[doublecount_mask]
    return i, j, D, S


def plot_atoms(atoms: Atoms,
               unitcell: np.ndarray | None = None) -> go.Figure:
    """Ball and stick plotly figure."""
    data = []

    # Atoms:
    points = SPHERE_POINTS
    triangles = triangulate_sphere()
    i, j, k = triangles.T
    for Z, xyz in zip(atoms.numbers, atoms.positions):
        x, y, z = (points * covalent_radii[Z] * 0.5 + xyz).T
        mesh = go.Mesh3d(x=x, y=y, z=z,
                         i=i, j=j, k=k, opacity=1,
                         color=color(Z))
        data.append(mesh)

    # Bonds:
    i, j, D, S = get_bonds(atoms)
    p = atoms.positions

    xyz = np.empty((3, len(i) * 3))
    xyz[:, 0::3] = p[i].T
    D[S.any(1)] *= 0.5
    xyz[:, 1::3] = (p[i] + D).T
    xyz[:, 2::3] = np.nan

    x, y, z = xyz
    data.append(go.Scatter3d(x=x, y=y, z=z, mode='lines',
                             line=dict(color='grey', width=10),
                             showlegend=False))

    # Unit cell:
    if unitcell is None:
        unitcell = atoms.cell
    x, y, z = (UNITCELL @ unitcell).T
    data.append(go.Scatter3d(x=x, y=y, z=z, mode='lines',
                             line=dict(color='#fa9fb5', width=6),
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
