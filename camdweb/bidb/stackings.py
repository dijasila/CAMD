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
            e = getattr(bilayer, 'binding_energy_gs', None)
            if e is None:
                e = bilayer.binding_energy_zscan
            rows.append((e, uid))
        tbl = table(
            ['Stacking',
             'Binding energy [meV/Å<sup>2</sup>]'],
            [[f'<a href="{uid}">{uid}</a>' if uid != material.uid else
              f'{uid}',
              f'{e:.3f}']
             for e, uid in sorted(rows)])
        pngfile = monolayer.folder / 'stackings.png'
        if not pngfile.is_file():
            create_figure(bilayers, pngfile)
        return PanelData(
            HTML.format(
                table=tbl,
                image=image(pngfile)),
            title='Stackings')


def create_figure(bilayers: dict[str, Material],
                  path: Path) -> None:
    fig, ax = plt.subplots()
    x = [bilayer.distance for bilayer in bilayers.values()]
    y = []
    for bilayer in bilayers.values():
        e = getattr(bilayer, 'binding_energy_gs', None)
        if e is None:
            e = bilayer.binding_energy_zscan
        y.append(e)
    ax.plot(x, y, 'o')
    ax.set_xlabel('distance [Å]')
    ax.set_ylabel('binding energy [meV/Å/Å]')
    plt.savefig(path)
    plt.close(fig)
