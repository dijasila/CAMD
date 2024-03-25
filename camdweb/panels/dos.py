from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData
from camdweb.html import image

HTML = """
<div class="row">
{img}
</div>
"""


class DOSPanel(Panel):
    def get_data(self,
                 material: Material) -> PanelData:
        return PanelData(
            HTML.format(
                img=image(material.folder / 'dos.png',
                          alt=f'DOS for {material.uid}')),
            title='Density of states')
