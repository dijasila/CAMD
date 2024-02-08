"""Panel base class."""
from __future__ import annotations
from typing import Callable, TYPE_CHECKING, Generator
if TYPE_CHECKING:
    from camdweb.material import Material


class Panel:
    title: str
    info = ''
    datafiles = []
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def __init__(self,
                 column_descriptions,
                 html_formatters):
        self.column_descriptions = column_descriptions
        self.html_formatters = html_formatters

    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        return

    def get_columns(self, material):
        return {}

    def table_rows(self, names, material):
        rows = []
        for name in names:
            value = material.columns.get(name)
            if value is not None:
                formatter = self.html_formatters.get(name, default_formatter)
                rows.append([self.column_descriptions.get(name, name),
                             formatter(value)])
        return rows


def default_formatter(value):
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        return f'{value:.3f}'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    return str(value)
