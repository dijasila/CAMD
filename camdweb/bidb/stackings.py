from pathlib import Path

import matplotlib.pyplot as plt
from camdweb.html import image, table
from camdweb.materials import Material
from camdweb.panels.panel import Panel, PanelData

HTML = """
<div class="row">
  <div class="col-6">
   {image}
  </div>
  <div class="col-6">
   {table}
  </div>
</div>
"""


class StackingsPanel(Panel):
    def get_data(self,
                 material: Material) -> PanelData:
        bilayers = material.data.get('bilayers')
        if bilayers is None:
            monolayer = material.data['monolayer']
            bilayers = monolayer.data['bilayers']
        else:
            monolayer = material
        rows = []
        for uid, bilayer in bilayers.items():
            e = bilayer.binding_energy_gs
            rows.append([f'<a href="{uid}">{uid}</a>',
                         f'{e:.3f}'])
        tbl = table(
            ['Stacking',
             'Binding energy [meV/Å<sup>2</sup>]'],
            rows)
        pngfile = monolayer.folder / 'stackings.png'
        if not pngfile.is_file():
            create_figure(bilayers, pngfile)
        return PanelData(
            HTML.format(
                table=tbl,
                image=image(pngfile)),
            title='Stackings')


def create_figure(bilayers: list[Material],
                  path: Path) -> None:
    fig, ax = plt.subplots()
    x = [bilayer.distance for bilayer in bilayers]
    y = [bilayer.binding_energy_gs for bilayer in bilayers]
    ax.plot(x, y)
    ax.set_xlabel('distance [Å]')
    ax.set_ylabel('binding energy [eV]')
    plt.savefig(path)
    fig.close()
