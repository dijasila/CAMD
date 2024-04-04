from typing import Generator

from ase.io import read
from asrlib.postprocessing.structureinfo import StructureInfo

from camdweb.panels.panel import Panel
from camdweb.material import Material
from camdweb.html import table


HTML = """
<div class="row">
{table}
</div>
"""


class SymmetryPanel(Panel):
    title = 'Crystal Symmetry'
    datafiles = ['structure.xyz']

    # SPGLIB Symmetry Analysis params
    symprec = 0.3
    angle_tolerance = 0.01
    hall_number = 0
    wrap = True

    # Cell Params
    # Layer Group Params
    # Formula Params

    def get_html(self, material: Material) -> Generator[str, None, None]:
        path = material.folder / self.datafiles[0]
        atoms = read(path)
        sinfo = StructureInfo(atoms=atoms)
        # generate spglib_data
        spglib_analysis = sinfo.get_spglib_analysis(
            symprec=self.symprec, angle_tolerance=self.angle_tolerance,
            hall_number=self.hall_number, wrap=self.wrap)
        # formula information
        formula = sinfo.formula(mode='metal', empirical=False)
        stoichimetry = sinfo.get_anonymous_reduced_formula
        crystal_type = spglib_analysis.crystal_type

        html = HTML.format(
            table=table(['#', 'Chemical symbol', 'Charges [|e|]'],
                        [(n, s, f'{c:.2f}') for n, (s, c)
                         in enumerate(zip(material.atoms.symbols,
                                          charges))]))
        # if (atoms.pbc == [True, True, False]).all():
        #     info['cell_area'] = abs(np.linalg.det(atoms.cell[:2, :2]))


