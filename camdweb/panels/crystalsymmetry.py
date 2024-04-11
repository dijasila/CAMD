from ase.io import read

from camdweb.panels.panel import Panel, PanelData
from camdweb.material import Material
from camdweb.html import table, image
from camdweb.oqmd12345.utils import StructureInfo

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


class CrystalSymmetryPanel(Panel):
    title = 'Crystal Symmetry'
    datafiles = ['structure.xyz']
    header_spgkey_lookup = {
        'pointgroup': 'Point group', 'wyckoffs': 'Wyckoff sites',
        'international': 'Space group', 'number': 'Space group No.',
        'layergroup': 'Layer group', 'lgnum': 'Layer group No.'}
    spglib_params = {'symprec': 0.3, 'angle_tolerance': 0.01, 'hall_number': 0,
                     'wrap': True}

    def get_data(self, material: Material) -> PanelData:
        col1 = self.create_column_one(material)
        col2 = self.create_column_two(material)
        return PanelData(
            HTML.format(col1=col1, col2=col2),
            title=self.title,
        )

    def create_column_one(self, material: Material) -> str:
        path = material.folder / self.datafiles[0]
        atoms = read(path)
        sinfo = StructureInfo(atoms=atoms)
        spglib_analysis = sinfo.get_spglib_analysis(**self.spglib_params)

        tableinfo = {header: spglib_analysis.dataset[spg_key]
                     for spg_key, header in self.header_spgkey_lookup.items()
                     if spg_key in spglib_analysis.dataset}
        return table(['Property', 'Value'],
                     list(map(lambda kvp: list(kvp), tableinfo.items())))

    def create_column_two(self, material: Material) -> str:
        col2 = image(material.folder / 'dos.png',  # change to BravisLattice
                     alt=f'Bravais Lattice for {material.uid}')
        return col2