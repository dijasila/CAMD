from typing import Callable
from cxdb.material import Material


class Panel:
    title: str
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def get_html(self, material: Material) -> tuple[str, str]:
        raise NotImplementedError

    def get_column_data(self, material: Material) -> dict[str, tuple[Any, str]]:
        raise NotImplementedError

    def update_column_data(self, material: Material) -> None:
        data = self.get_column_data(material)
        for name, (value, html) in data.items():
            material.add_column(name, value, html)
