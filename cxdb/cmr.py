import importlib
import sys
from pathlib import Path

from ase.db import connect
from bottle import Bottle

from cxdb.atoms import AtomsPanel
from cxdb.material import Material, Materials
from cxdb.utils import table
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
        return table(
            ['Project',
             'Number of materials',
             'Data',
             'Description'],
            [[project.title,
              len(project.materials),
              f'<a href="/{name}">{name}</a>',
              f'https://cmr.fysik.dtu.dk/{name}/{name}.html']
             for name, project in self.projects.items()])

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


class CMRAtomsPanel(AtomsPanel):
    def __init__(self, ndims: int):
        super().__init__(ndims)
        self.column_names.update(
            magstate='Magnetic state',
            ehull='Energy above convex hull [eV/atom]',
            hform='Heat of formation [eV/atom]',
            gap='Band gap (PBE) [eV]',
            energy='Energy [eV]',
            has_inversion_symmetry='Inversion symmetry',
            uid0='Old uid',
            evac='Vacuum level [eV]')
        self.columns = list(self.column_names)

    def update_data(self, material: Material):
        super().update_data(material)
        data = json.loads((material.folder / 'data.json').read_text())
        for key, value in data.items():
            material.add_column(key, value)


class ProjectDescription:
    def __init__(self, title, uid, ndims, column_names, initial_columns):
        self.title = title
        self.uid = uid
        self.ndims = ndims
        self.column_names = column_names
        self.initial_columns = initial_columns

    @classmethod
    def from_cmr_name(cls, name):
        mod = importlib.import_module(f'cmr.{name}.custom')
        Template = mod.Template
        column_names = {
            name: long + (f' [{unit}]' if unit else '')
            for name, (short, long, unit)
            in Template.raw_key_descriptions.items()}
        return cls(Template.title,
                   getattr(Template, 'uid', 'id'),
                   None,
                   column_names,
                   Template.default_columns)

    @classmethod
    def generic(cls, name):
        return cls(name, 'id', None, {}, ['formula'])


def app_from_db(db,
                project_description):
    pd = project_description
    root = Path()  # not used
    rows = []
    for row in db.select():
        atoms = row.toatoms()
        material = Material(root, str(row[pd.uid]), atoms)
        for name in pd.column_names:
            value = row.get(name)
            if value is not None:
                material.add_column(name, value)
        rows.append(material)

    ndims = sum(atoms.pbc) if pd.ndims is None else pd.ndims
    panels = [CMRAtomsPanel(ndims)]
    materials = Materials(rows, panels)
    return CMRProjectApp(materials, pd.initial_columns)


def main(filenames: list[str]) -> CMRProjectsApp:
    projects = {}
    for filename in filenames:
        path = Path(filename)
        name = path.stem
        project_description = ProjectDescription.from_cmr_name(name)
        db = connect(path)
        app = app_from_db(db, project_description)
        projects[name] = app

    return CMRProjectsApp(projects)


if __name__ == '__main__':
    main(sys.argv[1:]).app.run(host='0.0.0.0', port=8080, debug=True)
