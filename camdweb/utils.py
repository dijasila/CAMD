from __future__ import annotations
import re

import numpy as np
from ase.data import chemical_symbols


def formula_dict_to_string(count: dict[str, int]) -> str:
    """Convert dict to string representation.

    >>> formula_dict_to_strings({'H': 2, 'O': 1})
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
    H<sub>2</sub>O
    """
    return re.sub(r'(\d)', '<sub>\1</sub>', f)


def fft(atomic_numbers: list[int] | np.ndarray) -> tuple[dict[str, int],
                                                         str, str, str]:
    """Fast formula-transformations.

    Return dict mapping chemical symbols to number of chemical symbols.
    Three dicts are returned:

    >>> fft([1, 1, 8, 1, 1, 8])
    ({'O': 2, 'H': 4}, {'O': 1, 'H': 2}, {'A': 1, 'B': 2})

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


def doi(id: str | None) -> str | None:
    if id is None:
        return None
    assert isinstance(id, str)
    return f'<a href="https://doi.org/{id}">{id}</a>'


def cod(id: int | None) -> str | None:
    if id is None:
        return None
    assert isinstance(id, int)
    return f'<a href="{COD}/{id}.html">COD {id}</a>'


def icsd(id: int | None) -> str | None:
    if id is None:
        return None
    assert isinstance(id, int)
    return f'<a href="{ICSD}">ICSD {id}</a>'
