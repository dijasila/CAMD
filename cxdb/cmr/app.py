from __future__ import annotations

import argparse
import pickle
import sys
from functools import partial
from math import isfinite
from pathlib import Path
from typing import Callable

import rich.progress as progress
from ase.db import connect
from ase.db.row import AtomsRow
from bottle import Bottle, static_file, template

from cxdb.cmr.projects import ProjectDescription, create_project_description
from cxdb.html import FormPart, table
from cxdb.material import Material, Materials
from cxdb.panels.atoms import AtomsPanel
from cxdb.panels.panel import Panel
from cxdb.web import CXDBApp

CMR = 'https://cmr.fysik.dtu.dk'


class CMRProjectsApp:
    def __init__(self, project_apps: dict[str, CMRProjectApp]):
        self.project_apps = project_apps
        self.app = Bottle()
        self.app.route('/')(self.overview)
        self.app.route('/favicon.ico')(self.favicon)
        self.app.route('/<project_name>')(self.index1)
        self.app.route('/<project_name>/row/<uid>')(self.material)
        self.app.route('/<project_name>/callback')(self.callback)
        self.app.route('/<project_name>/png/<uid>')(self.png)

        for fmt in ['xyz', 'cif', 'json']:
            self.app.route(f'/<project_name>/<uid>/download/{fmt}')(
                partial(self.download, fmt=fmt))

    def download(self, project_name: str, uid: str, fmt: str):
        return self.project_apps[project_name].download(uid=uid, fmt=fmt)

    def overview(self) -> str:
        tbl = table(
            ['Project',
             'Number of materials',
             'Download data',
             'Description'],
            [[f'<a href="/{name}">{app.title}</a>',
              len(app.materials),
              (f'<a download="{name}.db" ' +
               f'href="{CMR}/_downloads/{name}.db">{name}.db</a>'),
              f'<a href="{CMR}/{name}/{name}.html">{name}</a>']
             for name, app in sorted(self.project_apps.items())])
        return template('cmr/overview', table=tbl, title='CMR projects')

    def index1(self, project_name: str) -> str:
        html = self.project_apps[project_name].index()
        return html.replace('/material/', f'/{project_name}/row/')

    def material(self, project_name: str, uid: str) -> str:
        html = self.project_apps[project_name].material(uid)
        html = html.replace('/callback', f'/{project_name}/callback')
        return html.replace('href="/">Search<',
                            f'href="/{project_name}">Search<', 1)

    def callback(self, project_name: str):
        return self.project_apps[project_name].callback()

    def favicon(self) -> bytes:
        path = Path(__file__).with_name('favicon.ico')
        return static_file(path.name, path.parent)

    def png(self, project_name: str, uid: str) -> bytes:
        app = self.project_apps[project_name]
        return app.png(uid, f'{project_name}/{uid}.png')


class CMRProjectApp(CXDBApp):
    def __init__(self,
                 materials: Materials,
                 initial_columns: list[str],
                 dbpath: Path,
                 title: str,
                 form_parts: list[FormPart]):
        super().__init__(materials, initial_columns)
        self.dbpath = dbpath
        self.title = title
        self.form_parts += form_parts

    def route(self) -> None:
        pass


class CMRAtomsPanel(AtomsPanel):
    def __init__(self,
                 column_names: dict[str, str],
                 create_column_one: Callable[[Material, Materials],
                                             tuple[str, str]],
                 create_column_two: Callable[[Material, Materials],
                                             tuple[str, str]]):
        super().__init__()
        self.column_names.update(column_names)
        self.column_names.update(
            energy='Energy [eV]',
            fmax='Maximum force [eV/Å]',
            smax='Maximum stress component [eV/Å<sup>3</sup>]',
            magmom='Total magnetic moment [μ<sub>B</sub>]')
        self.columns = list(self.column_names)
        self._create_column_one = create_column_one
        self._create_column_two = create_column_two

    def create_column_one(self, material, materials):
        col1, foot = self._create_column_one(material, materials)
        if not col1:
            return super().create_column_one(material, materials)
        return col1, foot

    def create_column_two(self, material, materials):
        col2, foot = self._create_column_two(material, materials)
        if not col2:
            return super().create_column_two(material, materials)
        return col2, foot


def app_from_db(dbpath: Path,
                project_description: ProjectDescription) -> CMRProjectApp:
    pd = project_description
    root = dbpath.parent
    pickle_file = dbpath.with_suffix('.pckl')
    if pickle_file.is_file():
        with open(pickle_file, 'rb') as fd:
            materials = pickle.load(fd)
    else:
        rows = []
        db = connect(dbpath)
        with progress.Progress() as pb:
            pid = pb.add_task('Reading rows:', total=len(db))
            for row in db.select():
                rows.append(row2material(row, pd, root))
                pb.advance(pid)

        panels: list[Panel] = [
            CMRAtomsPanel(pd.column_names,
                          pd.create_column_one,
                          pd.create_column_two)]
        panels += pd.panels
        materials = Materials(rows, panels)

    initial_columns = [name for name in pd.initial_columns
                       if name in materials.column_names]
    return CMRProjectApp(materials, initial_columns,
                         dbpath, pd.title, pd.form_parts)


def row2material(row: AtomsRow,
                 pd: ProjectDescription,
                 root: Path) -> Material:
    atoms = row.toatoms()
    if pd.pbc is not None:
        atoms.pbc = pd.pbc
    material = Material(root, str(row[pd.uid]), atoms)

    energy = row.get('energy')
    if energy is not None:
        material.add_column('energy', energy)
    forces = row.get('forces')
    if forces is not None:
        material.add_column('fmax', (forces**2).sum(axis=1).max()**0.5)
    stress = row.get('stress')
    if stress is not None:
        material.add_column('smax', abs(stress).max())
    magmom = row.get('magmom')
    if magmom is not None:
        material.add_column('magmom', magmom)

    for name in pd.column_names:
        value = row.get(name)
        if value is not None:
            if isinstance(value, int) and pd.column_names[name][-1] == ']':
                # Column description is something like
                # "Description ... [unit]".  This means we
                # have an integer with a unit!
                value = float(value)
            elif isinstance(value, float) and not isfinite(value):
                continue
            material.add_column(name, value)
    pd.postprocess(material)
    return material


def main(argv: list[str] | None = None) -> CMRProjectsApp:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'db_name', nargs='+',
        help='Filename of CMR-project SQLite- database.')
    parser.add_argument('--pickle', action='store_true')
    args = parser.parse_args(argv)
    project_apps = {}
    for filename in args.db_name:
        path = Path(filename)
        print(path)
        name = path.stem
        project_description = create_project_description(name)
        app = app_from_db(path, project_description)
        project_apps[name] = app
        if args.pickle:
            with open(path.with_suffix('.pckl'), 'wb') as fd:
                pickle.dump(app.materials, fd)

    return CMRProjectsApp(project_apps)


if __name__ == '__main__':
    app = main(sys.argv[1:])
    app.app.run(host='0.0.0.0', port=8082, debug=True)
