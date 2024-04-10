"""QPOD web-app.

This module has code to convert ~cmr/C2DB/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see camdweb.c2db.asr_panel module).
"""
from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from camdweb.materials import Material, Materials
from camdweb.html import table
from camdweb.c2db.asr_panel import ASRPanel, read_result_file
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.panels.charge_neutrality import ChargeNeutralitySuperpanel
from camdweb.panels.slater_janak import SlaterJanakPanel
from camdweb.panels.defect_symmetry import DefectSymmetryPanel
from camdweb.web import CAMDApp

HTML = """
<h4>{defect_title}</h4>
<div class="row">
  <div class="col-6">
    {column1}
  </div>
  <div class="col-6">
    {column2}
  </div>
</div>
"""

class QPODAtomsPanel(AtomsPanel):
    title = 'Summary'

    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        col1 = self.create_column_one(material)
        col2 = self.create_column_two(material)
        host = material.columns['host_name']
        defect = material.columns['defect_name']
        charge = material.columns['charge_state']
        yield HTML.format(column1=col1,
                          column2=col2,
                          defect_title=f'{defect} in {host} {charge}')

    def create_column_one(self,
                          material: Material) -> str:
        html1 = table(['Pristine crystal info', ''],
                      self.table_rows(material,
                                      ['host_name', 'host_crystal', 'host_spacegroup',
                                       'host_pointgroup', 'host_hof', 'host_gap_pbe',
                                       'host_gap_hse', 'host_uid'])) # add c2db link to host_uid
        html2 = table(['Defect properties', ''],
                      self.table_rows(material,
                                      ['r_nn']))
        html3 = table(['Basic properties', ''],
                      self.table_rows(material,
                                      ['magnetic']))
        return '\n'.join([html1, html2, html3])


def main(root: Path) -> CAMDApp:
    """Create QPOD app."""
    mlist: list[Material] = []
    files = list(root.glob('*/*/charge*'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(files))
        for f in files:
            data_file = f / 'data.json'
            if data_file.is_file():
                with open(data_file, 'r') as json_file:
                    data = json.load(json_file)
                    uid = data.get('uid')
                    material = Material.from_file(f / 'structure.xyz', uid)
                    material.columns.update(data)
                    mlist.append(material)
            else:
                print(f'''Warning: data.json file not found 
                        or not readable in {f}''')
            pb.advance(pid)

    panels: list[Panel] = [QPODAtomsPanel(),
                           ChargeNeutralitySuperpanel(),
                           SlaterJanakPanel(),
                           DefectSymmetryPanel()]
    # for name in ['bandstructure',
    #              'phonons',
    #              'bader']:
    #     panels.append(ASRPanel(name))
    #panels.append(ShiftCurrentPanel())

    materials = Materials(mlist, panels)

    initial_columns = ['host_name', 'defect_name', 'charge_state', 'formula', 'uid']

        # def __init__(self):
        # super().__init__()
        # self.column_names.update(
        #     magstate='Magnetic state',
        #     host_name='Host crystal',
        #     defect_name='Defect',
        #     charge_state='Charge',
        #     host_crystal='Host crystal type',
        #     host_uid='Host C2DB link',
        #     host_spacegroup='Host space group',
        #     host_pointgroup='Host point group',
        #     host_hof='Heat of formation of host material [eV/atom]',
        #     host_gap_pbe='Band gap of host material (PBE) [eV]',
        #     host_gap_hse='Band gap of host material (HSE) [eV]',
        #     energy='Energy [eV]',
        #     r_nn='Defect-defect distance [Ang]')
        # self.columns = list(self.column_names)

    materials.column_descriptions.update(
        magstate='Magnetic state',
        host_name='Host crystal',
        defect_name='Defect',
        charge_state='Charge',
        host_crystal='Host crystal type',
        host_uid='Host C2DB link',
        host_spacegroup='Host space group',
        host_pointgroup='Host point group',
        host_hof='Heat of formation of host material [eV/atom]',
        host_gap_pbe='Band gap of host material (PBE) [eV]',
        host_gap_hse='Band gap of host material (HSE) [eV]',
        energy='Energy [eV]',
        r_nn='Defect-defect distance [Ang]')

    return CAMDApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8083, debug=True)
