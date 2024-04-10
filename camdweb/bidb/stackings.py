from camdweb.panels.panel import Panel, PanelData, SkipPanel
from camdweb.materials import Material, Materials
from camdweb.html import table, image


class StackingsPanel(Panel):
    def get_data(self,
                 material: Material) -> PanelData:
        bilayers = material.data.get('bilayers')
        if bilayers is None:
            raise SkipPanel
        rows = []
        for uid, bilayer in bilayers.items():
            e = bilayer.binding_energy_gs
            rows.append([f'<a href="{uid}">{uid}</a>', f'{e:.3f}'])
        tbl = table(['Stacking', 'Binding energy [meV/Å<sup>2</sup>]'], rows)
        return PanelData(tbl, title='Stackings')


