"""Panel base class."""
from __future__ import annotations
import abc
from typing import Callable, TYPE_CHECKING, Generator
if TYPE_CHECKING:
    from camdweb.material import Material


class Panel(abc.ABC):
    title: str

    callbacks: dict[str, Callable[[Material, int], str]] = {}

    @abc.abstractmethod
    def get_html(self,
                 material: Material) -> Generator[str, None, None]:
        raise NotImplementedError

    def xxxupdate_data(self, material: Material) -> None:
        pass
