HTML = """
<img alt="DOS for {id}" src="/png/{id}/dos.png" />
"""


class DOSSection:
    def __init__(self):
        self.callbacks = {}

    def get_html(self, material):
        return ('Density of state',
                HTML.format(id=material.id),
                '')
