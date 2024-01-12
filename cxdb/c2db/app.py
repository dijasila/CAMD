"""C2DB web-app.

This module has code to convert ~cmr/C2DB/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see cxdb.asr_panel module).
"""
from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from cxdb.material import Material, Materials
from cxdb.panels.asr_panel import ASRPanel
from cxdb.panels.atoms import AtomsPanel
from cxdb.panels.panel import Panel
from cxdb.panels.shift_current import ShiftCurrentPanel
from cxdb.web import CXDBApp


class C2DBAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()
        self.column_names.update(
            magstate='Magnetic state',
            ehull='Energy above convex hull [eV/atom]',
            hform='Heat of formation [eV/atom]',
            gap='Band gap (PBE) [eV]',
            energy='Energy [eV]',
            has_inversion_symmetry='Inversion symmetry',
            uid0='Old uid',
            evac='Vacuum level [eV]')
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            material.add_column(key, value)


def main(root: Path) -> CXDBApp:
    """Create C2DB app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(files))
        for f in files:
            uid = f'{f.parent.name}-{f.name}'
            mlist.append(Material.from_file(f / 'structure.xyz', uid))
            pb.advance(pid)

    panels: list[Panel] = [C2DBAtomsPanel()]
    for name in ['bandstructure',
                 'phonons',
                 'bader']:
        panels.append(ASRPanel(name))
    panels.append(ShiftCurrentPanel())

    materials = Materials(mlist, panels)

    initial_columns = ['magstate', 'ehull', 'hform', 'gap', 'formula', 'area']

    return CXDBApp(materials, initial_columns, root)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8081, debug=True)
