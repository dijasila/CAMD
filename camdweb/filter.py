"""Code for efficient filtering of rows.

Example filter string:

* ``Cu>1``: more than 1 Cu atom
* ``gap>1.1``: "gap" larger than 1.1
* ``xc=PBE``: "xc" equal to "PBE"
* ``MoS2``: one or more MoS2 formula units

Strings can be combined with ``,`` (and) and ``|`` (or).  Grouping can be
done using ``(`` and ``)``:

* ``(Cu>1, gap>1.1) | Fe=0``
"""
from __future__ import annotations

import functools
import re
import sys
from collections import defaultdict
from typing import Callable

import numpy as np
from ase.data import chemical_symbols
from ase.formula import Formula

from camdweb import ColVal


@functools.lru_cache
def parse(q: str) -> Callable[[Index], set[int]]:
    """Convert filter string to Python function.

    >>> f = parse('H2,xc=PBE')
    >>> i = Index([('H', {'H': 2}, {'xc': 'PBE'})])
    Rows: 1 | Strings: 1 | Integers: 1 | Floats: 0 | Int-floats: 0
    >>> f(i)
    {0}
    >>> f = parse('gap > 5.0')
    >>> i = Index([('H', {'H': 2}, {'gap': 10.0}),
    ...            ('Si', {'Si': 2}, {'gap': 1.1})])
    Rows: 2 | Strings: 0 | Integers: 2 | Floats: 1 | Int-floats: 0
    >>> f(i)
    {0}
    """
    q = parse1(q)
    return eval(f'lambda i: {q}')


def parse1(q: str) -> str:
    """Quick'n'dirty hacky parsing of filter string to Python expression.

    In the Python expression, *i* is an Index object.

    >>> parse1('H2')
    "(i.formula('H2'))"
    >>> parse1('H2,xc=PBE')
    "(i.formula('H2')) & (i.key('xc', '=', 'PBE'))"
    >>> parse1('H2,Ag=0')
    "(i.formula('H2')) & (i.key('Ag', '=', 0))"
    >>> parse1('H>7,(Ag=1|Cu=1)')
    "(i.key('H', '>', 7)) & ((i.key('Ag', '=', 1)) | (i.key('Cu', '=', 1)))"
    """
    q0 = q
    q = q.replace(' ', '')
    if not q:
        return 'i.ids'
    n1 = len(q)
    p1 = q.count('(')
    p2 = q.count(')')
    if p1 != p2:
        raise SyntaxError('Mismatched parenthesis')

    # number of parsed characters (should match n1 when we are done):
    n = p1 + p2

    q = q.replace(',', ' & ')
    n2 = len(q)
    n += (n2 - n1) // 2

    q = q.replace('|', ' | ')
    n3 = len(q)
    n += (n3 - n2) // 2

    h: list[str] = []  # value strings
    for x in ['=', '<', '>', '<=', '>=', '!=']:
        # Find Ag=7, Cu<=7, ...
        matches = reversed(list(re.finditer(fr'([A-Z][a-z]?{x}[0-9]+)', q)))
        for m in matches:
            i, j = m.span()
            k, v = m[1].split(x)
            if k not in chemical_symbols:
                raise SyntaxError(f'Unknown chemical symbol "{k}"')
            q = q[:i] + f'(i.key({k!r}, {x!r}, {v}))' + q[j:]
            n += j - i

        # Find key=value, key>value, ...
        matches = reversed(
            list(re.finditer(fr'([A-Za-z][a-z0-9_]*{x}[-A-Za-z0-9.+_]+)', q)))
        for m in matches:
            i, j = m.span()
            k, v = m[1].split(x)
            if not k.islower():
                raise SyntaxError(f'Illegal key: "{k}"')
            v = repr(str2obj(v))
            q = q[:i] + f'(i.key({k!r}, {x!r}, #{len(h)}))' + q[j:]
            h.append(v)
            n += j - i

    # Find Ag7O14, ...
    matches = reversed(list(re.finditer(r'([A-Z][a-z]?[0-9]*)+', q)))
    for m in matches:
        i, j = m.span()
        if j < len(q) and q[j] == "'":
            continue  # we already dealt with this one
        q = q[:i] + f'(i.formula({m[0]!r}))' + q[j:]
        n += j - i

    # Insert value strings:
    for i, v in reversed(list(enumerate(h))):
        q = q.replace(f'#{i}', v)

    if n1 != n:
        raise SyntaxError(f'Bad filter string: {q0!r}')

    return q


def str2obj(s: str) -> bool | int | float | str:
    """Convert string to object.

    >>> str2obj('Hello')
    'Hello'
    >>> str2obj('7.4')
    7.4
    """
    x = {'True': True, 'False': False}.get(s)
    if x is not None:
        return x
    for type in [int, float]:
        try:
            return type(s)
        except ValueError:
            pass
    return s


class Index:
    def __init__(self,
                 rows: list[tuple[str,
                                  dict[str, int],
                                  dict[str, ColVal]]],
                 max_int_range: int = 350):
        """Indices for speeding up row filtering.

        The *rows* argument is a list of three element tuples:

        1) Reduced formula string with elements sorted by abundance followed by
           alphabetic sorting.  Examples: water: 'OH2', rocksalt: 'ClNa'.

        2) Atom count as a dict: {'H': 2, 'O': 1}.

        3) Key-value pairs: {'ok': True, 'gap': 1.1, ...}.

        >>> rows = [('BN', {'B': 2, 'N': 2}, {'gap': 3.6}),
        ...         ('C', {'C': 2}, {'gap': 0.0})]
        >>> i = Index(rows)
        Rows: 2 | Strings: 0 | Integers: 3 | Floats: 1 | Int-floats: 0
        >>> i.float_key('gap', '>', 0.0)
        {0}
        """

        integers = defaultdict(list)
        floats = defaultdict(list)
        self.strings: defaultdict[str, dict[str, set[int]]] = defaultdict(dict)
        self.natoms: dict[int, int] = {}
        self.reduced: defaultdict[str, set[int]] = defaultdict(set)
        self.ids = set()

        for i, (reduced, count, keys) in enumerate(rows):
            self.natoms[i] = sum(count.values())
            self.reduced[reduced].add(i)
            self.ids.add(i)
            for symbol, n in count.items():
                integers[symbol].append((n, i))
            for name, value in keys.items():
                assert name.islower(), name
                if isinstance(value, str):
                    dct = self.strings[name]
                    if value not in dct:
                        dct[value] = {i}
                    else:
                        dct[value].add(i)
                elif isinstance(value, float):
                    floats[name].append((value, i))
                elif isinstance(value, (int, bool)):
                    integers[name].append((int(value), i))
                else:
                    raise ValueError

        self.integers = {}
        self.floats = {}
        ni = 0  # number of ints converted to floats
        for name, idata in integers.items():
            idata.sort()
            ids = np.array([i for value, i in idata], dtype=np.int32)
            indices = [0]
            nmin = idata[0][0]
            nmax = idata[-1][0]
            if nmax - nmin > max_int_range:
                # Avoid too wide range of integer-index
                values = np.array([value
                                   for value, i in idata], dtype=np.int32)
                self.floats[name] = (values, ids)
                ni += 1
                continue
            m = nmin
            for j, (n, i) in enumerate(idata):
                if n > m:
                    indices += [j] * (n - m)
                    m = n
            indices.append(j + 1)
            self.integers[name] = (
                nmin,
                nmax,
                np.array(ids, dtype=np.int32),
                np.array(indices, dtype=np.int32))

        for name, fdata in floats.items():
            assert name not in self.integers, name
            fdata.sort()
            ids = np.array([i for value, i in fdata], dtype=np.int32)
            values = np.array([value for value, i in fdata])
            self.floats[name] = (values, ids)

        print(f'Rows: {len(rows)} | '
              f'Strings: {len(self.strings)} | '
              f'Integers: {len(self.integers)} | '
              f'Floats: {len(self.floats) - ni} | '
              f'Int-floats: {ni}')

    def key(self,
            name: str,
            op: str,
            value: bool | int | float | str) -> set[int]:
        if name in chemical_symbols:
            n = value  # number of atoms
            assert isinstance(n, int)
            if name not in self.integers and name not in self.floats:
                if op == '=' and n != 0:
                    return set()
                if op == '<' and n == 0:
                    return set()
                if op == '!=' and n == 0:
                    return set()
                if op == '>' and n == 0:
                    return set()
                if op == '>=' and n > 0:
                    return set()
                return self.ids

        if name in self.integers:
            assert isinstance(value, (int, bool)), (name, value)
            return self.integer_key(name, op, value)

        if name in self.floats:
            assert isinstance(value, (int, float))
            return self.float_key(name, op, value)

        if name in self.strings:
            value = str(value)
            if op == '=':
                if value in self.strings[name]:
                    return self.strings[name][value]
                return set()
            if op == '!=':
                result = set()
                for val, ids in self.strings[name].items():
                    if val != value:
                        result.update(ids)
                return result
            raise ValueError

        return set()

    def float_key(self, name: str, op: str, value: float) -> set[int]:
        values, ids = self.floats[name]
        if op == '!=':
            return self.ids - self.float_key(name, '=', value)
        if op == '<=':
            value = np.nextafter(value, value + 1)
            op = '<'
        if op == '>':
            value = np.nextafter(value, value + 1)
            op = '>='
        j1 = np.searchsorted(values, value)
        n = len(values)
        if op == '=':
            result = set()
            while j1 < n and values[j1] == value:
                result.add(ids[j1])
                j1 += 1
            return result
        if op == '<':
            return set(ids[:j1])
        if op == '>=':
            return set(ids[j1:])
        raise SyntaxError

    def integer_key(self, name: str, op: str, n: int) -> set[int]:
        nmin, nmax, ids, indices = self.integers[name]
        if op == '=':
            if nmin <= n <= nmax:
                d = n - nmin
                j1, j2 = indices[d:d + 2]
                return set(ids[j1:j2])
            return set()
        if op == '!=':
            return (self.integer_key(name, '<=', n - 1) |
                    self.integer_key(name, '>', n))
        if op == '<':
            return self.integer_key(name, '<=', n - 1)
        if op == '>=':
            return self.integer_key(name, '>', n - 1)
        if op == '<=':
            if n < nmin:
                return set()
            n = min(n, nmax)
            d = n - nmin
            j = indices[d + 1]
            return set(ids[:j])
        if op == '>':
            if n >= nmax:
                return set()
            if n < nmin:
                n = nmin - 1
            d = n - nmin
            j = indices[d + 1]
            return set(ids[j:])
        assert False, op

    def formula(self, f: str) -> set[int]:
        formula = Formula(f)
        if len(formula) == 1:
            return self.key(f, '>', 0)
        stoichiometry, reduced, n = formula.stoichiometry()
        ids = self.reduced[str(reduced)]
        if n == 1:
            return ids
        m = len(formula)
        return {id for id in ids if self.natoms[id] % m == 0}


if __name__ == '__main__':
    print(parse1(' '.join(sys.argv[1:])))
