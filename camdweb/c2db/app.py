"""C2DB web-app.

The camdweb.c2db.copy_files module has code to convert ~cmr/C2DB-ASR/tree/
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
# import multiprocessing as mp
from pathlib import Path

import rich.progress as progress

from camdweb.c2db.asr_panel import ASRPanel
from camdweb.html import Range, RangeX, Select, table
from camdweb.materials import Material, Materials
from camdweb.panels.atoms import AtomsPanel
from camdweb.panels.bader import BaderPanel
from camdweb.panels.bandstructure import BandStructurePanel
from camdweb.panels.convex_hull import ConvexHullPanel
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
                                      ['layergroup', 'lgnum', 'lable',
                                       'cod_id', 'icsd_id', 'doi', 'olduid']))
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

    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(folders))
        for f in folders:
            uid = f'{f.parent.name}-{f.name}'
            material = Material.from_file(f / 'structure.xyz', uid)
            data = json.loads((f / 'data.json').read_text())
            material.columns.update(data)
            mlist.append(material)
            pb.advance(pid)

    pool = None  # mp.Pool(maxtasksperchild=100)

    def asr_panel(name):
        return ASRPanel(name, pool)

    panels: list[Panel] = [
        C2DBAtomsPanel(),
        ConvexHullPanel(
            sources={'OQMD': ('Bulk crystals from OQMD123',
                              f'<a href={OQMD}/{{uid}}>{{formula:html}}</a>'),
                     'C2DB': ('Monolayers from C2DB',
                              '<a href={uid}>{formula:html}</a>')}),
        asr_panel('stiffness'),
        asr_panel('phonons'),
        asr_panel('deformationpotentials'),
        BandStructurePanel(),
        asr_panel('pdos'),
        asr_panel('effective_masses'),
        asr_panel('hse'),
        asr_panel('gw'),
        asr_panel('borncharges'),
        asr_panel('shg'),
        asr_panel('polarizability'),
        asr_panel('infraredpolarizability'),
        asr_panel('raman'),
        asr_panel('bse'),
        BaderPanel(),
        asr_panel('piezoelectrictensor'),
        ShiftCurrentPanel()]

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
        iscd_id='ICSD id of parent bulk structure',
        doi='Reported DOI',
        lgnum='Layer group number',
        label='Structure origin',
        gap='Band gap [eV]',
        gap_hse='Band gap (HSE06) [eV]',
        gap_gw='Band gap (G₀W₀) [eV]')

    materials.html_formatters.update(
        cod_id=cod,
        icsd_id=icsd,
        doi=doi,
        olduid=olduid)

    initial_columns = ['formula', 'ehull', 'hform', 'gap', 'magstate',
                       'layergroup']

    root = folders[0].parent.parent.parent
    app = CAMDApp(materials, initial_columns, root)
    app.form_parts += [
        Select('Dynamically stable', 'dyn_stab',
               ['', 'True', 'False'], ['', 'Yes', 'No']),
        Range('Energy above convex hull [eV/atom]', 'ehull', nonnegative=True),
        Select('Magnetic', 'magstate', ['', 'NM', 'FM'], ['', 'No', 'Yes']),
        RangeX('Band gap range [eV]', 'bg',
               ['gap', 'gap_hse', 'gap_gw'], ['PBE', 'HSE06', 'GW'])]
    return app


def test():  # pragma: no cover
    app = main(['AB2'])
    app.material_page('1MoS2-3')


def create_app():  # pragma: no cover
    """Create the WSGI app."""
    return main([str(path) for path in Path().glob('A*/')]).app


def check_all(pattern: str):  # pragma: no cover
    """Generate png-files."""
    c2db = main([str(path) for path in Path().glob(pattern)])
    for material in c2db.materials:
        print(material.uid)
        c2db.material_page(material.uid)


if __name__ == '__main__':
    main().app.run(host='0.0.0.0', port=8081, debug=True)
