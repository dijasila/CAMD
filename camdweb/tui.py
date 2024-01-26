from textual.app import App
from textual.widgets import Header, Footer, Input, Select, DataTable


class CAMD(App):  # pragma: no cover
    BINDINGS = [('d', 'toggle_dark', 'Toggle dark mode')]

    def compose(self):
        yield Header()
        yield Input(placeholder='filter: Au>2,H=1 or Cu,(H|O)')
        yield Select([(s, s) for s in ['AB2', 'AB', 'A']],
                     prompt='Stoichiometry')
        yield DataTable()
        yield Footer()

    def on_mount(self):
        table = self.query_one(DataTable)
        table.add_columns('A', 'B', 'C', 'D')
        table.add_rows(['asdf asdfasdf qweqwer xzcxcb'.split()] * 20)

    def on_input_submitted(self, event):
        print(event)


if __name__ == '__main__':  # pragma: no cover
    CAMD().run()
