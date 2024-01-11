from __future__ import annotations
import sys
from pathlib import Path

from ase.db import connect
from bottle import Bottle, static_file

from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.utils import table
from cxdb.web import CXDBApp
from cxdb.cmr.projects import create_project_description, ProjectDescription


class CMRProjectsApp:
    def __init__(self, project_apps: dict[str, CMRProjectApp]):
        self.project_apps = project_apps
        self.app = Bottle()
        self.app.route('/')(self.overview)
        self.app.route('/favicon.ico')(self.favicon)
        self.app.route('/<project_name>')(self.index1)
        self.app.route('/<project_name>/row/<uid>')(self.material)
        self.app.route('/<project_name>/callback')(self.callback)
        self.app.route('/<project_name>/download')(self.download_db_file)

    def overview(self) -> str:
        CMR = 'https://cmr.fysik.dtu.dk'
        return table(
            ['Project',
             'Number of materials',
             'Download data',
             'Description'],
            [[f'<a href="/{name}">{app.title}</a>',
              len(app.materials),
              f'<a download="" href="/{name}/download">{name}.db</a>',
              f'<a href="{CMR}/{name}/{name}.html">{name}</a>']
             for name, app in self.project_apps.items()])

    def index1(self, project_name: str) -> str:
        html = self.project_apps[project_name].index()
        return html.replace('/material/', f'/{project_name}/row/')

    def material(self, project_name: str, uid: str) -> str:
        return self.project_apps[project_name].material(uid)

    def callback(self, project_name: str, query: dict | None = None):
        return self.project_apps[project_name].callback(query)

    def download_db_file(self, project_name: str) -> bytes:
        path = self.project_apps[project_name].dbpath
        return static_file(path.name, path.parent)

    def favicon(self) -> bytes:
        path = Path(__file__).with_name('favicon.ico')
        return static_file(path.name, path.parent)


class CMRProjectApp(CXDBApp):
    def __init__(self,
                 materials: Materials,
                 initial_columns: list[str],
                 dbpath: Path,
                 title: str,
                 search: tuple[list[str], str]):
        super().__init__(materials, initial_columns)
        self.dbpath = dbpath
        self.title = title
        self.search_words, self.search_html = search

    def get_filter_string(self, query: dict) -> str:
        filter = super().get_filter_string(query)
        for word in self.search_words:
            q = query.get(word, '')
            if q:
                filter += f',{word}={q}'
        return filter.lstrip(',')

    def index(self, query: dict | None = None) -> str:
        html = super().index()
        if not self.search_html:
            return html
        return html.replace('<!-- EXTRA SELECT BLOCKS -->',
                            self.search_html, 1)

    def route(self) -> None:
        pass


class CMRAtomsPanel(AtomsPanel):
    def __init__(self, ndims: int, column_names: dict[str, str]):
        super().__init__(ndims)
        self.column_names.update(column_names)
        self.columns = list(self.column_names)


def app_from_db(dbpath: Path,
                project_description: ProjectDescription) -> CMRProjectApp:
    pd = project_description
    root = dbpath.parent  # not used
    rows = []
    for row in connect(dbpath).select():
        atoms = row.toatoms()
        material = Material(root, str(row[pd.uid]), atoms)
        for name in pd.column_names:
            value = row.get(name)
            if value is not None:
                material.add_column(name, value)
        rows.append(material)

    ndims = sum(atoms.pbc) if pd.ndims is None else pd.ndims
    panels = [CMRAtomsPanel(ndims, pd.column_names)]
    materials = Materials(rows, panels)
    initial_columns = [name for name in pd.initial_columns
                       if name in materials.column_names]
    return CMRProjectApp(materials, initial_columns,
                         dbpath, pd.title, pd.search)


def main(filenames: list[str]) -> CMRProjectsApp:
    project_apps = {}
    for filename in filenames:
        path = Path(filename)
        name = path.stem
        project_description = create_project_description(name)
        app = app_from_db(path, project_description)
        project_apps[name] = app

    return CMRProjectsApp(project_apps)


if __name__ == '__main__':
    main(sys.argv[1:]).app.run(host='0.0.0.0', port=8080, debug=True)
