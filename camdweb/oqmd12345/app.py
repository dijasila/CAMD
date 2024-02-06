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
from ase.io import read

from camdweb.material import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.panel import Panel
from camdweb.web import CAMDApp

COLUMN_NAMES = dict(
    etot='Total Energy [eV]',
    magstate='Magnetic state',
    fmax='Maximum Force [eV/A²]',
    smax='Maximum Stress',
    gap='Band gap (PBE) [eV]',
    gap_dir='Direct Band Gap (PBE) [eV]',
    oqmd_entry_id='OQMD Ref. ID')


class OQMD12345Material(Material):
    def __init__(self, folder: Path, uid: str):
        data = json.loads((folder / 'data.json').read_text())
        super().__init__(folder, uid, read(folder / 'structure.xyz'))
        self.etot: float = data['etot']
        self.magstate: str = data['magstate']
        self.fmax: float = data['fmax']
        self.smax: float = data['smax']
        self.gap: float = data['gap']
        self.gap_dir: float = data['gap_dir']
        self.oqmd_entry_id: int = data['oqmd_entry_id']

    def get_columns(self):
        columns = super().get_columns()
        columns.update(
            etot=self.etot,
            magstate=self.magstate,
            fmax=self.fmax,
            smax=self.smax,
            gap=self.gap,
            gap_dir=self.gap_dir,
            oqmd_entry_id=self.oqmd_entry_id)
        return columns


class OQMD12345AtomsPanel(AtomsPanel):
    def create_column_one(self,
                          material: OQMD12345Material) -> str:
        rows = []
        for key, name in COLUMN_NAMES.items():
            if key != 'oqmd_entry_id':
                value = getattr(self, key)
                rows.append([name, material.html_format_column(key, value)])
            else:
                rows.append([name, value
            ('Link to old C2DB', f'<a href={OLD}/{old}>{old}</a>')]

        return table([title, ''],
                     [[key, value]
                      for key, value in rows if value is not None]))
        return '\n'.join(tables)


def main(root: Path) -> CAMDApp:
    """Create CRYSP app."""
    mlist: list[Material] = []
    files = list(root.glob('A*/*/*/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading materials:', total=len(files))
        for f in files:
            uid = f'{f.parent.name}-{f.name}'
            mlist.append(Material(f, uid))
            pb.advance(pid)

    panels: list[Panel] = [OQMDAtomsPanel()]

    materials = Materials(mlist, panels, column_names)
    initial_columns = ['formula', 'gap', 'gap_dir',
                       'magstate', 'etot', 'volume']
    return CAMDApp(materials, initial_columns, root)


def create_app():  # pragma: no cover
    return main(Path()).app


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8086, debug=True)
