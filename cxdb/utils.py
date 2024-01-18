import numpy as np
from ase.data import chemical_symbols


def formula_dict_to_strings(count: dict[str, int]) -> tuple[str, str]:
    """Convert dict to string representations: text and HTML.

    >>> formula_dict_to_strings({'H': 2, 'O': 1})
    ('H2O', 'H<sub>2</sub>O')
    """
    s = ''
    html = ''
    for symbol, c in count.items():
        s += symbol
        html += symbol
        if c > 1:
            s += str(c)
            html += f'<sub>{c}</sub>'
    return s, html


def fft(atomic_numbers: list[int] | np.ndarray) -> tuple[dict[str, int],
                                                         dict[str, int],
                                                         dict[str, int]]:
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
    return count, reduced, stoichiometry
