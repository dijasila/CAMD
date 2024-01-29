"""CRYSP web-app.

This module has code to convert ~crysp*/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see cxdb.asr_panel module).
"""
from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from camdweb.panels.atoms import AtomsPanel
from camdweb.material import Material, Materials
from camdweb.panels.panel import Panel
from camdweb.web import CAMDApp


class CRYSPAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()
        self.column_names.update(
            etot='Total Energy [eV]',
            magstate='Magnetic state',
            fmax='Maximum Force [eV/A²]',
            smax='Maximum Stress',
            gap='Band gap (PBE) [eV]',
            gap_dir='Direct Band Gap (PBE) [eV]',
            oqmd_entry_id='OQMD Ref. ID')
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            material.add_column(key, value)


def main(root: Path) -> CAMDApp:
    """Create CRYSP app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading materials:', total=len(files))
        for f in files:
            uid = f'{f.parent.name}-{f.name}'
            mlist.append(Material.from_file(f / 'structure.xyz', uid))
            pb.advance(pid)

    panels: list[Panel] = [CRYSPAtomsPanel()]
    # for name in ['bandstructure',
    #              'phonons',
    #              'bader']:
    #     panels.append(ASRPanel(name))

    # breakpoint()
    materials = Materials(mlist, panels)
    initial_columns = ['gap', 'formula']

    return CAMDApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8081, debug=True, server='waitress')
