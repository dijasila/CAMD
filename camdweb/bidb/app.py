from __future__ import annotations

import json
from pathlib import Path

import rich.progress as progress

from camdweb.html import table, image, Select, Range
from camdweb.materials import Material, Materials
from camdweb.panels.misc import MiscPanel
from camdweb.panels.atoms import AtomsPanel
from camdweb.web import CAMDApp
from camdweb.bidb.stackings import StackingsPanel
from camdweb.c2db.bs_dos_bz_panel import BSDOSBZPanel
from camdweb.c2db.asr_panel import ASRPanel, read_result_file

COLUMN_DESCRIPTIONS = {
    'binding_energy_zscan': 'Binding energy (zscan)',
    'number_of_layers': 'Number of layers',
    'monolayer_uid': 'Monolayer ID',
    'bilayer_uid': 'Bilayer ID',
    'dynamically_stable': 'Dynamically stable',
    'magnetic': 'Magnetic',
    'interlayer_magnetic_state': 'Interlayer magnetic state',
    'interlayer_magnetic_exchange': 'Interlayer magnetic exchange',
    'slide_stability': 'Slide Stability',
    'binding_energy_gs': 'Binding Energy (gs) [meV/Å<sup>2</sup>]',
    'ehull': 'Energy above convex hull [eV/atom]',
    'gap_pbe': 'Band gap (PBE)',
    'icsd_id': 'ICSD id of parent bulk structure',
    'cod_id': 'COD id of parent bulk structure',
    'layer_group': 'Layer group',
    'layer_group_number': 'Layer group number',
    'space_group': 'Space group',
    'space_group_number': 'Space group number',
    'distance': 'Distance [Å]'}


class BiDBAtomsPanel(AtomsPanel):
    column_descriptions = COLUMN_DESCRIPTIONS

    def update_material(self, material):
        data = json.loads((material.folder / 'data.json').read_text())
        material.columns.update(data)

    def create_column_one(self,
                          material: Material) -> str:
        skip = {'interlayer_magnetic_exchange',
                'binding_energy_zscan',
                'monolayer_uid'}
        keys = COLUMN_DESCRIPTIONS.keys() - skip
        rows = self.table_rows(material, keys)
        if material.number_of_layers == 2:
            uid = material.data['monolayer'].uid
            rows.append(
                ['Monolayer ID',
                 f'<a href="/material/{uid}">{material.monolayer_uid}</a>'])
        return table(None, rows)


def read_material(path: Path,
                  uid: str) -> Material:
    material = Material.from_file(path / 'structure.xyz', uid)
    try:
        evac = read_result_file(path / 'results-asr.gs.json')['evac']
    except FileNotFoundError:
        pass
    else:
        material.columns['evac'] = evac
    return material


def main(root: Path, pattern: str = '*') -> CAMDApp:
    mlist: list[Material] = []
    paths = list(root.glob(f'{pattern}/*/monolayer/'))
    with progress.Progress() as pb:
        pid = pb.add_task('Reading matrerials:', total=len(paths))
        for f1 in paths:
            uid1 = f1.parent.name
            monolayer = read_material(f1, uid1)
            bilayers = {}
            for f2 in f1.parent.glob('*/'):
                if f2.name != 'monolayer':
                    uid2 = f'{f2.parent.name}-{f2.name}'
                    bilayer = read_material(f2, uid2)
                    bilayer.data['monolayer'] = monolayer
                    bilayers[uid2] = bilayer
                    mlist.append(bilayer)
            monolayer.data['bilayers'] = bilayers
            mlist.append(monolayer)
            pb.advance(pid)

    panels = [BiDBAtomsPanel(),
              StackingsPanel(),
              BSDOSBZPanel(),
              ASRPanel('fermisurface'),
              ASRPanel('raman'),
              MiscPanel()]

    materials = Materials(mlist, panels)

    initial_columns = [
        'formula',
        'number_of_layers',
        'binding_energy_gs',
        'slide_stability',
        'uid',
        'magnetic']

    app = CAMDApp(materials, initial_columns, root=root)

    app.title = 'BiDB'
    app.logo = image('bidb-logo.png', alt='BiDB-logo')
    app.links = [
        ('CMR', 'https://cmr.fysik.dtu.dk'),
        ('BiDB', 'https://cmr.fysik.dtu.dk/bidb/bidb.html')]
    app.form_parts += [
        Select('Number of layers', 'number_of_layers', ['', '1', '2']),
        Range('Binding energy [meV/Å<sup>2</sup>] (bilayers)',
              'binding_energy_gs'),
        Select('Slide stability (bilayers)', 'slide_stability',
               ['', 'Stable']),
        Range('Band gap range [eV]', 'gap_pbe'),
        Select('Magnetic', 'magnetic', ['', '0', '1'])]

    return app


def create_app():
    """Create the WSGI app."""
    app = main(Path())
    return app.app


def check_all(pattern: str):  # pragma: no cover
    """Generate png-files."""
    bidb = main(Path(), pattern)
    for material in bidb.materials:
        print(material.uid)
        bidb.material_page(material.uid)


if __name__ == '__main__':
    main(Path()).app.run(host='0.0.0.0', port=8084, debug=True)
