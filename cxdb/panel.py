"""Panel base class."""
from __future__ import annotations
import abc
from typing import Callable, TYPE_CHECKING
if TYPE_CHECKING:
    from cxdb.material import Material, Materials


class Panel(abc.ABC):
    title: str
    column_names: dict[str, str] = {}

    callbacks: dict[str, Callable[[Material, int], str]] = {}

    @abc.abstractmethod
    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        raise NotImplementedError

    def update_data(self, material: Material) -> None:
        pass
