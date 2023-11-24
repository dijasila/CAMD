from typing import Callable
from cxdb.material import Material


class Section:
    title: str
    callbacks: dict[str, Callable[[Material, int], str]] = {}

    def get_html(self, material: Material) -> tuple[str, str]:
        raise NotImplementedError
