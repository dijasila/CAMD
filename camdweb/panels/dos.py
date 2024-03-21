from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData

HTML = """
<div class="row">
<img alt="DOS for {uid}" src="/png/{uid}/dos.png" />
</div>
"""


class DOSPanel(Panel):
    def get_data(self,
                 material: Material) -> PanelData:
        return PanelData(HTML.format(uid=material.uid),
                         title='Density of states')
