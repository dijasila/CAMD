from camdweb.material import Material
from camdweb.panels.panel import Panel, PanelData
from camdweb.html import table

HTML = """
<div class="row">
  <div class="col-6">
   {col1}
  </div>
  <div class="col-6">
   {col2}
  </div>
</div>
"""


class MiscPanel(Panel):
    def __init__(self):
      self.subpanels = None

    def get_data(self,
                 material: Material) -> PanelData:

        keys = []
        dkeys = material.columns.keys()
        for key in dkeys:
            keys.append(key)
        kl2 = len(keys) // 2
        kl1 = len(keys) - kl2

        html1 = table(['Miscellaneous details', ''],
                      self.table_rows(material,
                                      keys[0:kl1]))
        html2 = table(['Miscellaneous details', ''],
                      self.table_rows(material,
                                      keys[kl1:]))

        pd = PanelData(HTML.format(col1=html1, col2=html2),
                         title='Miscellaneous')

        if self.subpanels is not None:
            for subpanel in self.subpanels:
                pd.subpanels.append(subpanel.get_data(material))

        return pd
