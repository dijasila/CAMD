import sys
from pathlib import Path

from ase.db import connect
from bottle import Bottle, static_file

from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.utils import table
from cxdb.web import CXDBApp
from cxdb.cmr_projects import create_project_description


class CMRProjectsApp:
    def __init__(self, projects):
        self.projects = projects
        self.app = Bottle()
        self.app.route('/')(self.overview)
        self.app.route('/<project_name>')(self.index)
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
            [[f'<a href="/{name}">{project.title}</a>',
              len(project.materials),
              f'<a diwnload="" href="/{name}/download">{name}.db</a>',
              f'<a href="{CMR}/{name}/{name}.html">{name}</a>']
             for name, project in self.projects.items()])

    def index(self, project_name) -> str:
        html = self.projects[project_name].index()
        return html.replace('/material/', f'/{project_name}/row/')

    def material(self, project_name, uid) -> str:
        return self.projects[project_name].material(uid)

    def callback(self, project_name):
        return self.projects[project_name].callback()

    def download_db_file(self, project_name: str) -> bytes:
        path = self.projects[project_name].dbpath
        return static_file(path.name, path.parent)


class CMRProjectApp(CXDBApp):
    def __init__(self, materials, dbpath, title, initial_columns):
        super().__init__(materials, initial_columns)
        self.dbpath = dbpath
        self.title = title

    def route(self):
        pass


class CMRAtomsPanel(AtomsPanel):
    def __init__(self, ndims: int, column_names: dict[str, str]):
        super().__init__(ndims)
        self.column_names.update(column_names)
        self.columns = list(self.column_names)


def app_from_db(dbpath,
                project_description):
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
    return CMRProjectApp(materials, dbpath, pd.title, initial_columns)


def main(filenames: list[str]) -> CMRProjectsApp:
    projects = {}
    for filename in filenames:
        path = Path(filename)
        name = path.stem
        project_description = create_project_description(name)
        app = app_from_db(path, project_description)
        projects[name] = app

    return CMRProjectsApp(projects)


if __name__ == '__main__':
    main(sys.argv[1:]).app.run(host='0.0.0.0', port=8080, debug=True)
