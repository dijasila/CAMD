"""Panel base class."""
from __future__ import annotations
from typing import Callable, TYPE_CHECKING, Generator
from camdweb import ColVal
if TYPE_CHECKING:
    from camdweb.material import Material


class Panel:
    title: str
    info = ''
    datafiles: list[str] = []
    column_descriptions: dict[str, str] = {}
    html_formatters: dict[str, Callable[[ColVal, bool], str]] = {}
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        yield 'hello'

    def get_columns(self, material):
        return {}

    def update_material(self, material):
        pass

    def update_column_descriptions(self, column_descriptions):
        column_descriptions.update(self.column_descriptions)
        self.column_descriptions = column_descriptions

    def update_html_formatters(self, html_formatters):
        html_formatters.update(self.html_formatters)
        self.html_formatters = html_formatters

    def table_rows(self, names, material):
        rows = []
        for name in names:
            value = material.columns.get(name)
            if value is not None:
                formatter = self.html_formatters.get(name, default_formatter)
                rows.append([self.column_descriptions.get(name, name),
                             formatter(value, link=True)])
        return rows


def default_formatter(value, link=False):
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        return f'{value:.3f}'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    return str(value)
