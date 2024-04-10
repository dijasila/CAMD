"""Panel base class."""
from __future__ import annotations

import abc
<<<<<<< HEAD
from typing import TYPE_CHECKING, Callable, Generator, Iterable, List

=======
from typing import TYPE_CHECKING, Callable, Iterable
from dataclasses import dataclass
>>>>>>> origin/main
from camdweb import ColVal

if TYPE_CHECKING:
    from camdweb.material import Material


class SkipPanel(Exception):
    """Don't show this panel."""


@dataclass
class PanelData:
    html: str
    title: str
    info: str = ''
    script: str = ''
    subpanels: list[PanelData] | None = None


class Panel(abc.ABC):
    title: str = "Unnamed"
    info = ''
    datafiles: list[str] = []
    column_descriptions: dict[str, str] = {}
    html_formatters: dict[str, Callable[..., str]] = {}
    callbacks: dict[str, Callable[[Material, int], str]] = {}
<<<<<<< HEAD
    
    def __init__(self):
        super().__init__()
        self.subpanels = None
=======

    def __init__(self):
        super().__init__()
>>>>>>> origin/main

    @abc.abstractmethod
    def get_data(self,
                 material: Material) -> PanelData:
        raise NotImplementedError

    def update_material(self, material: Material) -> None:
        pass

    def update_column_descriptions(self,
                                   column_descriptions: dict[str, str]
                                   ) -> None:
        column_descriptions.update(self.column_descriptions)
        self.column_descriptions = column_descriptions

    def update_html_formatters(self,
                               html_formatters: dict[
                                   str,
                                   Callable[..., str]]
                               ) -> None:
        html_formatters.update(self.html_formatters)
        self.html_formatters = html_formatters

    def table_rows(self,
                   material: Material,
                   names: Iterable[str]) -> list[list[str]]:
        rows = []
        for name in names:
            value = material.columns.get(name)
            if value is not None:
                formatter = self.html_formatters.get(name, default_formatter)
                rows.append([self.column_descriptions.get(name, name),
                             formatter(value, link=True)])
        return rows

<<<<<<< HEAD
    def add_subpanels(self, material: Material): 
        return None

    def generate_webpanel(self, material: Material):
        self.add_subpanels(material)               
        self.generator = self.get_html(material)

        if self.subpanels is not None:
            for subpanel in self.subpanels:
                subpanel.generate_webpanel(material = material)
        try:
            html = next(self.generator)
        except StopIteration:
            self.generator = None
            return
        if not html == '':  # result is ready
            self.generator = iter([html])

    def get_webpanel(self): # To be called after generation started.
        assert self.generator is not None
        
        html = next(self.generator)
        html, script = cut_out_script(html)
        subwebpanels = None

        if self.subpanels is not None:
            subwebpanels = list()
            for subpanel in self.subpanels:
                subwebpanels.append(subpanel.get_webpanel())

        # Should i purge the generators here?
        self.generator = None

        return WebPanel(self.title, self.info, html, subwebpanels), script

def cut_out_script(html: str) -> tuple[str, str]:
    r"""We need to put the script tags in the footer.

    >>> cut_out_script('''Hello
    ... <script>
    ...   ...
    ... </script>
    ... CAMd''')
    ('Hello\n\nCAMd', '<script>\n  ...\n</script>')
    """
    m = re.search(r'(<script.*</script>)', html, re.MULTILINE | re.DOTALL)
    if m:
        i, j = m.span()
        return html[:i] + html[j:], html[i:j]
    return html, ''            
=======
>>>>>>> origin/main

def default_formatter(value: ColVal, link: bool = False) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        return f'{value:.3f}'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    return str(value)
<<<<<<< HEAD

# Simple class for parsing panel and subpanel html info to the front end.
class WebPanel:
    def __init__(self, panel_title, info, html, subpanels = None):
        self.panel_title = panel_title
        self.info = info
        self.html = html
        self.subpanels = subpanels

    def get_properties(self):
        subpanel_properties = None
        if self.subpanels is not None:
            subpanel_properties = list()
            for subpanel in self.subpanels:
                subpanel_properties.append(subpanel.get_properties())
        return self.panel_title, self.info, self.html, subpanel_properties
=======
>>>>>>> origin/main
