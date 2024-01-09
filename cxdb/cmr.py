import sys
from pathlib import Path

from ase.db import connect
from bottle import Bottle

from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.web import CXDBApp


class CMRProjectsApp:
    def __init__(self, projects):
        self.projects = projects
        self.app = Bottle()
        self.app.route('/')(self.overview)
        self.app.route('/<project_name>')(self.index)
        self.app.route('/<project_name>/row/<uid>')(self.material)
        self.app.route('/<project_name>/callback')(self.callback)

    def overview(self):
        return '\n'.join(self.projects)

    def index(self, project_name):
        html = self.projects[project_name].index()
        return html.replace('/material/', f'/{project_name}/row/')

    def material(self, project_name, uid):
        return self.projects[project_name].material(uid)

    def callback(self, project_name):
        return self.projects[project_name].callback()


class CMRProjectApp(CXDBApp):
    def __init__(self, materials, initial_columns):
        super().__init__(materials, initial_columns)

    def route(self):
        pass


def app_from_db(db,
                uid='id',
                initial_columns=None,
                ndims: int | None = None):
    root = Path()  # not used
    rows = []
    for row in db.select():
        atoms = row.toatoms()
        if ndims is None:
            ndims = sum(atoms.pbc)
        material = Material(root, str(row[uid]), atoms)
        rows.append(material)
    panels = [AtomsPanel(ndims)]
    materials = Materials(rows, panels)
    initial_columns = initial_columns or ['formula']
    return CMRProjectApp(materials, initial_columns)


def main(filenames: list[str]) -> CMRProjectsApp:
    projects = {}
    for filename in filenames:
        path = Path(filename)
        name = path.stem
        db = connect(path)
        app = app_from_db(db)
        projects[name] = app

    return CMRProjectsApp(projects)


if __name__ == '__main__':
    main(sys.argv[1:]).app.run(host='0.0.0.0', port=8080, debug=True)
