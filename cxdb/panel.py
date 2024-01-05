import abc
from typing import Callable
from cxdb.material import Material, Materials


class Panel(abc.ABC):
    title: str
    column_names: dict[str, str] = {}

    callbacks: dict[str, Callable[[Material, int], str]] = {}

    @abc.abstractmethod
    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        ...

    def update_data(self, material: Material) -> None:
        pass
