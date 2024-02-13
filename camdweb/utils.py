from __future__ import annotations

import multiprocessing as mp
import re
from contextlib import contextmanager

import numpy as np
from ase.data import chemical_symbols
from camdweb import ColVal


def formula_dict_to_string(count: dict[str, int]) -> str:
    """Convert dict to string representation.

    >>> formula_dict_to_string({'H': 2, 'O': 1})
    'H2O'
    """
    s = ''
    for symbol, c in count.items():
        s += symbol
        if c > 1:
            s += str(c)
    return s


def html_format_formula(f: str) -> str:
    """Convert formula string to HTML.

    >>> html_format_formula('H2O')
    'H<sub>2</sub>O'
    """
    return re.sub(r'(\d+)', r'<sub>\1</sub>', f)


def fft(atomic_numbers: list[int] | np.ndarray) -> tuple[dict[str, int],
                                                         str, str, str]:
    """Fast formula-transformations.

    Return dict mapping chemical symbols to number of chemical symbols.
    Three dicts are returned:

    >>> fft([1, 1, 8, 1, 1, 8])
    ({'O': 2, 'H': 4}, 'O2H4', 'OH2', 'AB2')

    full, reduced and stoichiometry.

    (We need to do this a million times and the ase.formula module is
    too slow).
    """
    values, counts = np.unique(atomic_numbers, return_counts=True)
    nunits = np.gcd.reduce(counts)
    symbols = [chemical_symbols[v] for v in values]
    count = {symbol: c
             for c, symbol in sorted((c, symbol)
                                     for symbol, c in zip(symbols, counts))}
    reduced = {symbol: c // nunits for symbol, c in count.items()}
    abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    stoichiometry = {abc[i]: c
                     for i, (symbol, c) in enumerate(reduced.items())}
    return (count,
            formula_dict_to_string(count),
            formula_dict_to_string(reduced),
            formula_dict_to_string(stoichiometry))


COD = 'https://www.crystallography.net/cod/'
ICSD = 'https://icsd.products.fiz-karlsruhe.de/en/'


def doi(id: ColVal, link: bool = False) -> str:
    """Create HTML anchor link to DOI.

    >>> doi('10.1103/ABC.95.216', True)
    '<a href="https://doi.org/10.1103/ABC.95.216">10.1103/ABC.95.216</a>'
    >>> doi('10.1103/ABC.95.216', False)
    '10.1103/ABC.95.216'
    """
    assert isinstance(id, str)
    if link:
        return f'<a href="https://doi.org/{id}">{id}</a>'
    return id


def cod(id: ColVal, link: bool = False) -> str:
    assert isinstance(id, (str, int))
    if link:
        return f'<a href="{COD}/{id}.html">COD {id}</a>'
    return str(id)


def icsd(id: ColVal, link: bool = False) -> str:
    assert isinstance(id, (str, int))
    if link:
        return f'<a href="{ICSD}">ICSD {id}</a>'
    return str(id)


class NoPool:
    def imap_unordered(self, worker, work):
        for args in work:
            yield worker(args)


@contextmanager
def process_pool(processes: int):
    if processes == 1:
        yield NoPool()
    else:
        yield mp.Pool(processes=processes)
