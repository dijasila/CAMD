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
from cxdb.panel import Panel
from cxdb.utils import table

HTML = """
<h4>Basic properties: {formula}</h4>
<div class="row">
  <div class="col-6">
   {table}
  </div>
  <div class="col-6">
    <label>Repeat:</label>
    <select onchange="cb(this.value, 'atoms', '{uid}')">
      <option value="1">1</option>
      <option value="2">2</option>
      <option value="3">3</option>
    </select>
    <div id='atoms' class='atoms'></div>
    {axes}
  </div>
</div>
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

    def __init__(self, ndims: int) -> None:
        self.callbacks = {'atoms': self.plot}
        self.ndims = ndims
        self.column_names = {}
        self.columns = ['stoichiometry', 'uid']
        if ndims:
            if ndims == 1:
                unit = 'Å'
            else:
                unit = f'Å<sup>{ndims}</sup>'
            name = DIMS[ndims]
            self.column_names[name] = f'{name.title()} [{unit}]'
            self.columns.append(name)

    def update_data(self, material):
        if self.ndims > 0:
            pbc = material.atoms.pbc
            vol = abs(np.linalg.det(material.atoms.cell[pbc][:, pbc]))
            material.add_column(DIMS[self.ndims], vol)

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        tbl = table(None,
                    [(materials.column_names[name], material[name])
                     for name in self.columns
                     if name in material._values])
        return (HTML.format(table=tbl,
                            axes=self.axes(material),
                            uid=material.uid,
                            formula=material['formula']),
                FOOTER.format(atoms_json=self.plot(material, 3)))

    def axes(self, material: Material) -> str:
        atoms = material.atoms
        tbl1 = table(
            ['Axis', 'x [Å]', 'y [Å]', 'y [Å]', 'Periodic'],
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


UNITCELL = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1],
            [nan, nan, nan], [1, 0, 0], [1, 0, 1],
            [nan, nan, nan], [0, 1, 0], [0, 1, 1],
            [nan, nan, nan], [1, 1, 0], [1, 1, 1]]


def plot_atoms(atoms: Atoms,
               unitcell: np.ndarray | None = None) -> go.Figure:
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
    p = atoms.positions
    i, j, D, S = neighbor_list('ijDS', atoms,
                               cutoff=covalent_radii[atoms.numbers] * 1.2)
    xyz = np.empty((3, len(i) * 3))
    xyz[:, 0::3] = p[i].T
    D[S.any(1)] *= 0.5
    xyz[:, 1::3] = (p[i] + D).T
    xyz[:, 2::3] = np.nan
    x, y, z = xyz
    data.append(go.Scatter3d(x=x, y=y, z=z, mode='lines'))

    # Unit cell:
    if unitcell is None:
        unitcell = atoms.cell
    x, y, z = (UNITCELL @ unitcell).T
    data.append(go.Scatter3d(x=x, y=y, z=z, mode='lines'))

    fig = go.Figure(data=data)
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.update_xaxes(showgrid=False)
    fig.update_layout(template='simple_white')
    fig.update_scenes(aspectmode='data')
    return fig


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
    return hull.simplices.copy()


@cache
def color(Z: int) -> str:
    a, b, c = (jmol_colors[Z] * 255).astype(int)
    return f'rgb({a},{b},{c})'


if __name__ == '__main__':
    atoms = read(sys.argv[1])
    assert isinstance(atoms, Atoms)
    fig = plot_atoms(atoms)
    fig.show()
