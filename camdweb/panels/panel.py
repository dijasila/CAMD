"""Panel base class."""
from __future__ import annotations

import abc
from typing import TYPE_CHECKING, Callable, Iterable
from dataclasses import dataclass, field
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
    subpanels: list[PanelData] = field(default_factory=list)

    def build_sidebar(self, super_pid: int=None):
        from bottle import template

        if len(self.subpanels) == 0 and super_pid is None:
            panel_case = 'Panel'
        elif super_pid is None:
            panel_case = 'SuperPanel'
        elif super_pid is not None and len(self.subpanels) == 0:
            panel_case = 'SubPanel'
        else:
            raise ValueError('Sub-panels do not support sub-panels.')

        sidebar_html = template(
            'sidebar_element.html',
            case=panel_case,
            title=self.title,
            pid=id(self),
            ppid=super_pid)

        for subpanel in self.subpanels:
            subpanel_html = subpanel.build_sidebar(super_pid=id(self))
            sidebar_html = '\n'.join([sidebar_html, subpanel_html])

        return sidebar_html

    def build_body(self, super_pid: int=None):
        from bottle import template

        if len(self.subpanels) == 0 and super_pid is None:
            panel_case = 'Panel'
        elif super_pid is None:
            panel_case = 'SuperPanel'
        elif super_pid is not None and len(self.subpanels) == 0:
            panel_case = 'SubPanel'
        else:
            raise ValueError('Sub-panels do not support sub-panels.')

        body_html = template(
            'panel_element.html',
            case=panel_case,
            title=self.title,
            html=self.html,
            info=self.info,
            pid=id(self),
            ppid=super_pid)

        for subpanel in self.subpanels:
            subpanel_html = subpanel.build_body(super_pid=id(self))
            body_html = '\n'.join([body_html, subpanel_html])

        return body_html


class Panel(abc.ABC):
    title: str = "Unnamed"
    info = ''
    datafiles: list[str] = []
    column_descriptions: dict[str, str] = {}
    html_formatters: dict[str, Callable[..., str]] = {}
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def __init__(self):
        super().__init__()

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


def default_formatter(value: ColVal, link: bool = False) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        return f'{value:.3f}'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    return str(value)
