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
from ase.io import read
from camdweb.c2db.asr_panel import ASRPanel
from camdweb.html import Range, RangeX, Select, table
from camdweb.material import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.bader import BaderPanel
from camdweb.panels.bandstructure import BandStructurePanel
from camdweb.panels.convex_hull import ConvexHullPanel
from camdweb.panels.panel import Panel
from camdweb.panels.shift_current import ShiftCurrentPanel
from camdweb.utils import cod, doi, icsd
from camdweb.web import CAMDApp

OLD = 'https://cmrdb.fysik.dtu.dk/c2db/row/'


class C2DBMaterial(Material):
    def __init__(self, folder: Path, uid: str):
        super().__init__(folder, uid, read(folder / 'structure.xyz'))
        data = json.loads((folder / 'data.json').read_text())
        self.olduid: str = data['olduid']
        self.has_inversion_symmetry: bool = data['has_inversion_symmetry']
        self.gap: float = data['gap']
        self.evac: float = data['evac']
        self.hform: float = data['hform']
        self.magstate: str = data['magstate']
        self.ehull: float = data['ehull']
        self.energy: float = data['energy']
        self.spin_axis: str = data['spin_axis']
        self.efermi: float = data['efermi']
        self.dyn_stab: bool = data['dyn_stab']
        self.layergroup: str = data['layergroup']
        self.lgnum: int = data['lgnum']
        self.gap_hse: float | None = data.get('gap_pbe')
        self.gap_gw: float | None = data.get('gap_gw')
        self.cod_id: int | None = data.get('cod_id')
        self.icsd_id: int | None = data.get('icsd_id')
        self.doi: str | None = data.get('doi')
        self.label: str | None = data.get('label')


class C2DBAtomsPanel(AtomsPanel):
    def __init__(self):
        super().__init__()

    def create_column_one(self,
                          material: C2DBMaterial,
                          materials: Materials) -> str:
        old = material.olduid
        table1 = [
            ('Layer group', material.layergroup),
            ('Layer group number', material.lgnum),
            ('COD id of parent bulk structure', cod(material.cod_id)),
            ('ICSD id of parent bulk structure', icsd(material.icsd_id)),
            ('Structure origin', material.label),
            ('Reported DOI', doi(material.doi)),
            ('Link to old C2DB', f'<a href={OLD}/{old}>{old}</a>')]

        table2 = [
            ('Energy above convex hull [eV/atom]', material.ehull),
            ('Heat of formation [eV/atom]', material.hform),
            ('Dynamically stable', 'Yes' if material.dyn_stab else 'No')]

        table3 = [
            ('Magnetic', 'Yes' if material.magstate == 'FM' else 'No'),
            ('Band gap [eV]', material.gap),
            ('Band gap (HSE06) [eV]', material.gap_hse),
            ('Band gap (G₀W₀) [eV]', material.gap_gw)]

        # return three tables with None values removed:
        tables = []
        for title, rows in [('Structure info', table1),
                            ('Stability', table2),
                            ('Basic properties', table3)]:
            tables.append(
                table([title, ''],
                      [[key, value]
                       for key, value in rows if value is not None]))
        return '\n'.join(tables)


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
            material = C2DBMaterial(f, uid)
            mlist.append(material)
            pb.advance(pid)

    panels: list[Panel] = [
        C2DBAtomsPanel(),
        ConvexHullPanel(),
        ASRPanel('stiffness', keys),
        ASRPanel('phonons', keys),
        ASRPanel('deformationpotentials', keys),
        ASRPanel('bandstructure', keys),
        BandStructurePanel(),
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

    column_names = dict(
        has_inversion_symmetry='Inversion symmetry',
        gap='Band gap (PBE) [eV]',
        evac='Vacuum level [eV]',
        hform='Heat of formation [eV/atom]',
        olduid='Old uid',
        magstate='Magnetic state',
        ehull='Energy above convex hull [eV/atom]',
        energy='Energy [eV]',
        spin_axis='Spin axis',
        efermi='Fermi level [eV]',
        dyn_stab='Dynamically stable')

    materials = Materials(mlist, panels, column_names)

    initial_columns = ['formula', 'ehull', 'hform', 'gap', 'magstate', 'area']

    root = folders[0].parent.parent.parent
    app = CAMDApp(materials, initial_columns, root)
    app.form_parts += [
        Select('Dynamically stable', 'dyn_stab',
               ['', 'True', 'False'], ['', 'Yes', 'No']),
        Range('Energy above convex hull [eV/atom]', 'ehull'),
        Select('Magnetic', 'magstate', ['', 'NM', 'FM'], ['', 'No', 'Yes']),
        RangeX('Band gap range [eV]', 'bg',
               ['gap', 'gap_hse', 'gap_gw'], ['PBE', 'HSE06', 'GW'])]
    return app


def test():  # pragma: no cover
    app = main(['AB2'])
    app.material_page('1MoS2-2')


def create_app():  # pragma: no cover
    return main([path.name for path in Path().glob('A*/')]).app


if __name__ == '__main__':
    main().app.run(host='0.0.0.0', port=8081, debug=True)
