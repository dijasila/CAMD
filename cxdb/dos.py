from cxdb.panel import Panel
from cxdb.material import Material, Materials

HTML = """
<div class="row">
<img alt="DOS for {uid}" src="/png/{uid}/dos.png" />
</div>
"""


class DOSPanel(Panel):
    title = 'Density of states'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        return (HTML.format(uid=material.uid), '')
