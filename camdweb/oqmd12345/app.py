"""OQMD12345 web-app.

This module has code to convert ~crysp*/tree/ folders and friends
(see PATTERNS variable below) to canonical tree layout.

::

    python -m camdweb.oqmd12345.app
"""
from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from camdweb.html import table
from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.web import CAMDApp

COLUMN_DESCRIPTIONS = dict(
    etot='Total Energy [eV]',
    magstate='Magnetic state',
    fmax='Maximum Force [eV/A²]',
    smax='Maximum Stress',
    gap='Band gap (PBE) [eV]',
    gap_dir='Direct Band Gap (PBE) [eV]',
    oqmd_entry_id='OQMD Ref. ID')

OQMD = 'https://oqmd.org/materials/entry'


class OQMD12345AtomsPanel(AtomsPanel):
    column_descriptions = COLUMN_DESCRIPTIONS.copy()

    def create_column_one(self,
                          material: Material) -> str:
        rows = self.table_rows(material, COLUMN_DESCRIPTIONS)
        return table(None, rows)


def oqmd(id, link=False):
    if link:
        return f'<a href={OQMD}/{id}>{id}</a>'
    return id


def main(root: Path) -> CAMDApp:
    """Create CRYSP app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/structure.xyz'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading materials:', total=len(files))
        for f in files:
            f = f.parent
            uid = f'{f.parent.name}-{f.name}'
            material = Material.from_file(f / 'structure.xyz', uid)
            data = json.loads((f / 'data.json').read_text())
            material.columns.update(data)
            mlist.append(material)
            pb.advance(pid)

    panels: list[Panel] = [OQMD12345AtomsPanel()]

    materials = Materials(mlist, panels)

    materials.html_formatters['oqmd_entry_id'] = oqmd

    initial_columns = ['formula', 'gap', 'gap_dir',
                       'magstate', 'etot', 'volume']
    return CAMDApp(materials, initial_columns, root)


def create_app():  # pragma: no cover
    return main(Path()).app


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8086, debug=True)
