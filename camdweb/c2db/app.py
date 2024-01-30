"""C2DB web-app.

The camdweb.c2db.copy_files module has code to convert ~cmr/C2DB-ASR/tree/
folders and friends (see PATTERNS variable below) to a canonical tree
layout.

Also contains simple web-app that can run off the tree of folders.

The goal is to have the code decoupled from ASE, GPAW, CMR and ASR.
Right now ASR webpanel() functions are still used
(see camdweb.c2db.asr_panel module).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import rich.progress as progress
from camdweb.material import Material, Materials
from camdweb.c2db.asr_panel import ASRPanel
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.panels.bader import BaderPanel
from camdweb.panels.shift_current import ShiftCurrentPanel
from camdweb.panels.convex_hull import ConvexHullPanel
from camdweb.web import CAMDApp


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


def main(argv: list[str] | None = None) -> CAMDApp:
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
                if key == 'uid':
                    assert value == uid
                    continue  # already added
                material.add_column(key, value)
                keys.add(key)
            pb.advance(pid)

    panels: list[Panel] = [
        C2DBAtomsPanel(),
        ConvexHullPanel(),
        ASRPanel('stiffness', keys),
        ASRPanel('phonons', keys),
        ASRPanel('deformationpotentials', keys),
        ASRPanel('bandstructure', keys),
        # ASRPanel('pdos', keys),
        ASRPanel('effective_masses', keys),
        ASRPanel('hse', keys),
        ASRPanel('gw', keys),
        ASRPanel('borncharges', keys),
        ASRPanel('shg', keys),
        ASRPanel('polarizability', keys),
        ASRPanel('infraredpolarizability', keys),
        ASRPanel('raman', keys),
        # ASRPanel('bse', keys),
        BaderPanel(),
        ASRPanel('piezoelectrictensor', keys),
        ShiftCurrentPanel()]

    materials = Materials(mlist, panels)

    initial_columns = ['formula', 'ehull', 'hform', 'gap', 'magstate', 'area']

    root = folders[0].parent.parent.parent
    return CAMDApp(materials, initial_columns, root)


def test():  # pragma: no cover
    app = main(['AB2'])
    app.material('1MoS2-1')


def create_app():
    return main([path.name for path in Path().glob('A*/')]).app


if __name__ == '__main__':
    main().app.run(host='0.0.0.0', port=8081, debug=True)
