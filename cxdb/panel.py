from typing import Callable, Any
from cxdb.material import Material


def creates(*filenames):
    def deco(func1):
        def func2(material):
            if (material.folder / filenames[0]).is_file():
                return
            func1(material)
        return func2
    return deco


class Panel:
    title: str
    column_names: dict[str, str] = {}
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def get_html(self,
                 material: Material,
                 column_names: dict[str, str]) -> tuple[str, str]:
        raise NotImplementedError

    def get_column_data(self,
                        material: Material) -> dict[str, tuple[Any, str]]:
        return {}

    def update_column_data(self, material: Material) -> None:
        data = self.get_column_data(material)
        for name, (value, html) in data.items():
            material.add_column(name, value, html)
