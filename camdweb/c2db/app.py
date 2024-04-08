"""C2DB web-app.

The camdweb.c2db.copy module has code to convert ~cmr/C2DB-ASR/tree/
folders and friends (see PATTERNS variable below) to a canonical tree
layout.

Also contains web-app that can run off the tree of folders.

The goal is to have the code decoupled from ASE, GPAW, CMR and ASR.
Right now ASR webpanel() functions are still used
(see camdweb.c2db.asr_panel module).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import rich.progress as progress

from camdweb.c2db.asr_panel import ASRPanel
from camdweb.c2db.bs_dos_bz_panel import BSDOSBZPanel
from camdweb.html import Range, RangeX, Select, table, image
from camdweb.materials import Material, Materials
from camdweb.optimade.app import add_optimade
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.bader import BaderPanel
from camdweb.panels.convex_hull import ConvexHullPanel
from camdweb.panels.emass import EmassPanel
from camdweb.panels.panel import Panel
from camdweb.panels.shift_current import ShiftCurrentPanel
from camdweb.utils import cod, doi, icsd
from camdweb.web import CAMDApp

OLD = 'https://cmrdb.fysik.dtu.dk/c2db/row/'
OQMD = 'https://cmrdb.fysik.dtu.dk/oqmd123/row'


class C2DBAtomsPanel(AtomsPanel):
    def create_column_one(self,
                          material: Material) -> str:
        html1 = table(['Structure info', ''],
                      self.table_rows(material,
                                      ['layergroup', 'lgnum', 'label',
                                       'cod_id', 'icsd_id', 'doi']))
        html2 = table(['Stability', ''],
                      self.table_rows(material,
                                      ['ehull', 'hform', 'dyn_stab']))
        html3 = table(['Basic properties', ''],
                      self.table_rows(material,
                                      ['magstate', 'gap', 'gap_hse',
                                       'gap_gw']))
        return '\n'.join([html1, html2, html3])


def olduid(uid, link=False):  # pragma: no cover
    if link:
        return f'<a href={OLD}/{uid}>{uid}</a>'
    return uid


class C2DBApp(CAMDApp):
    """C2DB app with /row/<olduid> endpoint."""

    title = 'C2DB'
    logo = image('c2db-logo.png', alt='C2DB-logo')
    links = [
        ('CMR', 'https://cmr.fysik.dtu.dk'),
        ('C2DB', 'https://cmr.fysik.dtu.dk/c2db/c2db.html')]

    def __init__(self,
                 materials: Materials,
                 initial_columns: list[str],
                 root: Path | None = None,
                 olduid2uid: dict[str, str] | None = None):
        super().__init__(materials,
                         initial_columns=initial_columns,
                         initial_filter_string='dyn_stab=True, ehull<0.2',
                         root=root)
        self.olduid2uid = olduid2uid or {}

    def route(self):
        super().route()
        self.app.route('/row/<olduid>')(self.material_page_old_uid)

    def material_page_old_uid(self, olduid: str) -> str:
        return self.material_page(self.olduid2uid[olduid])


def main(argv: list[str] | None = None) -> CAMDApp:
    """Create C2DB app."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'path', nargs='+',
        help='Path to tree.  Examples: "AB2" (same as "AB2/*/*"), '
        '"AB2/1MoS2" (same as "AB2/1MoS2/*/") or "AB2/1MoS2/1".')
    args = parser.parse_args(argv)

    folders = []
    for path in args.path:
        p = Path(path)
        if p.name.startswith('A'):
            folders += list(p.glob('*/*/'))
        elif p.name.rstrip('t').isdigit():  # 't' for temporary (AB2/1MoS2/1t)
            folders.append(p)
        else:
            folders += list(p.glob('*/'))

    olduid2uid = {}
    mlist: list[Material] = []
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            uid = f'{f.parent.name}-{f.name}'
            material = Material.from_file(f / 'structure.xyz', uid)
            data = json.loads((f / 'data.json').read_text())
            material.columns.update(data)
            mlist.append(material)
            if hasattr(material, 'olduid'):
                olduid2uid[material.olduid] = uid
            pb.advance(pid)

    panels: list[Panel] = [
        C2DBAtomsPanel(),
        ConvexHullPanel(
            sources={'OQMD': ('Bulk crystals from OQMD123',
                              f'<a href={OQMD}/{{uid}}>{{formula:html}}</a>'),
                     'C2DB': ('Monolayers from C2DB',
                              '{formula:html}, <a href={uid}>{uid}</a>')}),
        ASRPanel('stiffness'),
        ASRPanel('phonons'),
        ASRPanel('deformationpotentials'),
        BSDOSBZPanel(),
        EmassPanel(),
        ASRPanel('fermisurface'),
        ASRPanel('hse'),
        ASRPanel('gw'),
        ASRPanel('borncharges'),
        ASRPanel('shg'),
        Polarizability(),
        ASRPanel('infraredpolarizability'),
        ASRPanel('raman'),
        ASRPanel('bse'),
        BaderPanel(),
        ASRPanel('piezoelectrictensor'),
        ShiftCurrentPanel(),
        ASRPanel('collect_spiral'),
        ASRPanel('dmi'),
        ASRPanel('spinorbit')]

    materials = Materials(mlist, panels)

    materials.column_descriptions.update(
        has_inversion_symmetry='Inversion symmetry',
        evac='Vacuum level [eV]',
        hform='Heat of formation [eV/atom]',
        olduid='Old uid',
        magstate='Magnetic state',
        ehull='Energy above convex hull [eV/atom]',
        energy='Energy [eV]',
        spin_axis='Spin axis',
        efermi='Fermi level [eV]',
        dyn_stab='Dynamically stable',
        cod_id='COD id of parent bulk structure',
        icsd_id='ICSD id of parent bulk structure',
        doi='Mono/few-layer report(s)',
        layergroup='Layer group',
        lgnum='Layer group number',
        label='Structure origin',
        gap='Band gap [eV]',
        gap_hse='Band gap (HSE06) [eV]',
        gap_gw='Band gap (G₀W₀) [eV]',
        folder='Original file-system folder')

    materials.html_formatters.update(
        cod_id=cod,
        icsd_id=icsd,
        doi=doi,
        olduid=olduid)

    initial_columns = ['formula', 'ehull', 'hform', 'gap', 'magstate',
                       'layergroup']

    root = folders[0].parent.parent.parent
    app = C2DBApp(materials,
                  initial_columns,
                  root,
                  olduid2uid)
    app.form_parts += [
        Select('Magnetic:', 'magstate', ['', 'NM', 'FM'], ['-', 'No', 'Yes']),
        Select('Dynamically stable:', 'dyn_stab',
               ['', 'True', 'False'], ['-', 'Yes', 'No'], default='True'),
        Range('Energy above convex hull [eV/atom]:', 'ehull',
              nonnegative=True, default=('0.0', '0.2')),
        RangeX('Band gap range [eV]:', 'bg',
               ['gap', 'gap_hse', 'gap_gw'], ['PBE', 'HSE06', 'GW'])]

    return app


def test():  # pragma: no cover
    app = main(['AB2'])
    app.material_page('1MoS2-3')


def create_app():  # pragma: no cover
    """Create the WSGI app."""
    app = main([str(path) for path in Path().glob('A*/')])
    add_optimade(app)
    return app.app


def check_all(pattern: str):  # pragma: no cover
    """Generate png-files."""
    c2db = main([str(path) for path in Path().glob(pattern)])
    for material in c2db.materials:
        print(material.uid)
        c2db.material_page(material.uid)


def check_def_pot():  # pragma: no cover
    from .asr_panel import read_result_file
    for path in Path().glob('A*/*/*/results-asr.deformationpotentials.json'):
        r = read_result_file(path)
        try:
            r['defpots_soc']
        except KeyError:
            f = json.loads(path.with_name('data.json').read_text())['folder']
            print(f)


if __name__ == '__main__':
    main().app.run(host='0.0.0.0', port=8081, debug=True)
