from typing import Callable
from cxdb.material import Material, Materials


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
                 materials: Materials) -> tuple[str, str]:
        raise NotImplementedError

    def update_data(self, material: Material) -> None:
        pass
