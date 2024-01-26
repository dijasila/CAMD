from typing import Generator

from camdweb.material import Material, Materials
from camdweb.panels.panel import Panel

HTML = """
<div class="row">
<img alt="DOS for {uid}" src="/png/{uid}/dos.png" />
</div>
"""


class DOSPanel(Panel):
    title = 'Density of states'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> Generator[str, None, None]:
        yield HTML.format(uid=material.uid)
