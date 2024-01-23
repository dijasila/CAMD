"""C2DB web-app.

This module has code to convert ~cmr/C2DB/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

Also contains simple web-app that can run off the tree of folders.

Goal is to have the code decoupled from ASE, GPAW and ASR.
Right now ASR webpanel() functions are still used (see cxdb.asr_panel module).
"""
from __future__ import annotations

import argparse
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
            has_inversion_symmetry='Inversion symmetry',
            gap='Band gap (PBE) [eV]',
            evac='Vacuum level [eV]',
            hform='Heat of formation [eV/atom]',
            uid0='Old uid',
            magstate='Magnetic state',
            ehull='Energy above convex hull [eV/atom]',
            energy='Energy [eV]')
        self.columns = list(self.column_names)


def main(argv: list[str] | None = None) -> CXDBApp:
    """Create C2DB app."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'path', nargs='+',
        help='Path to tree.  Examples: "AB2" (same as "AB2/*/*"), '
        '"AB2/1MoS2" (same as "AB2/1MoS2/*/") or "AB2/1MoS2/1".')
    args = parser.parse_args(argv)

    mlist: list[Material] = []
    folders = []
    for path in args.path:
        p = Path(path)
        if p.name.startswith('A'):
            folders += list(p.glob('*/*/'))
        elif p.name.isdigit():
            folders.append(p)
        else:
            folders += list(p.glob('*/'))

    keys = set()
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            uid = f'{f.parent.name}-{f.name}'
            material = Material.from_file(f / 'structure.xyz', uid)
            mlist.append(material)
            data = json.loads((f / 'data.json').read_text())
            for key, value in data.items():
                material.add_column(key, value)
                keys.add(key)
            pb.advance(pid)

    panels: list[Panel] = [C2DBAtomsPanel()]
    for name in ['stiffness',
                 'phonons',
                 'deformationpotentials',
                 'bandstructure',
                 'pdos',
                 'effective_masses',
                 'hse',
                 'gw',
                 'borncharges',
                 'shg',
                 'polarizability',
                 'infraredpolarizability',
                 'raman',
                 # 'bse',
                 'bader',
                 'piezoelectrictensor']:
        print(name)
        panels.append(ASRPanel(name, keys))
    panels.append(ShiftCurrentPanel())

    materials = Materials(mlist, panels)

    initial_columns = ['formula', 'ehull', 'hform', 'gap', 'magstate', 'area']

    root = folders[0].parent.parent.parent
    return CXDBApp(materials, initial_columns, root)


def test():  # pragma: no cover
    app = main(['AB2'])
    app.material('2MoS2-2')


if __name__ == '__main__':
    main().app.run(host='0.0.0.0', port=8081, debug=True)
