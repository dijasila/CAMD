from cxdb.panel import Panel
from cxdb.material import Material

HTML = """
<img alt="DOS for {id}" src="/png/{id}/dos.png" />
"""


class DOSPanel(Panel):
    title = 'Density of states'

    def get_html(self, material: Material) -> tuple[str, str]:
        return (HTML.format(id=material.id), '')
