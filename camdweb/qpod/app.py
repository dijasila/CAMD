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

from camdweb.material import Material, Materials
from camdweb.c2db.asr_panel import ASRPanel, read_result_file
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.panels.charge_neutrality import ChargeNeutralityPanel
from camdweb.panels.slater_janak import SlaterJanakPanel
from camdweb.web import CAMDApp


class QPODAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()
        self.column_names.update(
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
            r_nn='Defect-defect distance [Ang]')  # change to r_nn
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            if key != 'uid':
                material.add_column(key, value)


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
                    mlist.append(Material.from_file(f / 'structure.xyz', 
                                                    uid=uid))
            else:
                print(f'''Warning: data.json file not found 
                        or not readable in {f}''')
            pb.advance(pid)

    panels: list[Panel] = [QPODAtomsPanel(),
                           ChargeNeutralityPanel(),
                           SlaterJanakPanel()]
    # for name in ['bandstructure',
    #              'phonons',
    #              'bader']:
    #     panels.append(ASRPanel(name))
    #panels.append(ShiftCurrentPanel())

    materials = Materials(mlist, panels)

    initial_columns = ['host_name', 'defect_name', 'charge_state', 'formula', 'uid']

    return CAMDApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8083, debug=True)
