from bottle import TEMPLATE_PATH, Bottle, request, static_file, template

from cxdb.web import CXDBApp


class CMRProjectsApp:
    def __init__(self, projects):
        self.projects = projects
        self.app = Bottle()
        self.app.route('/')(self.overview)
        self.app.route('/callback')(self.callback)
        self.app.route('/<project>')(self.index)
        self.app.route('/<project>/row/<uid>')(self.material)

    def material(self, project_name, uid):
        return self.projects[project_name].material(uid)


class CMRProjectApp(CXDBApp):
    def __init__(self):
        super().__init__()

    def route(self):
        pass

