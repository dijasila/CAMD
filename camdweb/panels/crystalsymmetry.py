from typing import Generator

from ase.io import read

from camdweb.panels.panel import Panel, PanelData
from camdweb.material import Material
from camdweb.html import table, image

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

    # SPGLIB Symmetry Analysis params
    symprec = 0.3
    angle_tolerance = 0.01
    hall_number = 0
    wrap = True

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
        from asrlib.postprocessing.structureinfo import StructureInfo
        sinfo = StructureInfo(atoms=atoms)

        # generate spglib_data
        spglib_analysis = sinfo.get_spglib_analysis(
            symprec=self.symprec, angle_tolerance=self.angle_tolerance,
            hall_number=self.hall_number, wrap=self.wrap)

        sinfo_keys = ['pointgroup', 'wyckoffs', 'international', 'number',
                      'layergroup', 'lgnum']
        header = ['Point group', 'Wyckoff sites', 'Space group',
                  'Space group No.', 'Layer group', 'Layer group No.']
        data = {key: spglib_analysis.dataset[dataset_key]
                for key, dataset_key in zip(header, sinfo_keys)
                if dataset_key in spglib_analysis.dataset}
        return table(['Property', 'Value'],
                     list(map(lambda kvp: list(kvp), data.items())))

    def create_column_two(self, material: Material) -> str:
        col2 = image(material.folder / 'dos.png',  # change to BravisLattice
                     alt=f'Bravais Lattice for {material.uid}')
        return col2