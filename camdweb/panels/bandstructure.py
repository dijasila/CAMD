from camdweb.panels.panel import Panel


class BandStructurePanel(Panel):
    title = 'Band structure'

    def get_html(self, material, materials):
        yield '¡¡HELLO BAND STRUCTURE!!'
