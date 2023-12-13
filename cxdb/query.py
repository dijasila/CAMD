from __future__ import annotations
import functools
import re
import sys
from collections import defaultdict
from typing import Callable

import numpy as np
from ase.data import chemical_symbols
from ase.formula import Formula


@functools.lru_cache
def parse(q: str) -> Callable[[Index], set[int]]:
    """Convert query string to Python function.

    >>> f = parse('H2,xc=PBE')
    >>> i = Index([({'H': 2}, {'xc': 'PBE'})])
    >>> f(i)
    {0}
    """
    q = parse1(q)
    return eval(f'lambda i: {q}')


def parse1(q: str) -> str:
    """Quick'n'dirty hacky parsing of query string to Python expression.

    In the Python expression, *n* is a dict mapping chemical symbols their
    numbers and *k* is a key-value dict:

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
            q = q[:i] + f'(i.key({k!r}, {x!r}, {v}))' + q[j:]
            n += j - i

        # Find key=value, key>value, ...
        matches = reversed(
            list(re.finditer(fr'([a-z][a-z0-9_]+{x}[-A-Za-z0-9.+_]+)', q)))
        for m in matches:
            i, j = m.span()
            k, v = m[1].split(x)
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
        raise SyntaxError(f'Bad query string: {q0!r}')

    return q


def str2obj(s: str) -> bool | int | float | str:
    """Convert string to object.

    >>> str2obj('Hello')
    'Hello'
    >>> str2obj('7.4')
    7.4
    """
    x = {'True': True, 'False': False}.get(s)
    if x:
        return x
    for type in [int, float]:
        try:
            return type(s)
        except ValueError:
            pass
    return s


class Index:
    def __init__(self, rows):
        integers = defaultdict(list)
        floats = defaultdict(list)
        self.strings = defaultdict(lambda: defaultdict(set))
        self.ids = set()
        for i, (count, keys) in enumerate(rows):
            self.ids.add(i)
            for symbol, n in count.items():
                integers[symbol].append((n, i))
            for name, value in keys.items():
                if isinstance(value, str):
                    self.strings[name][value].add(i)
                elif isinstance(value, float):
                    floats[name].append((value, i))
                elif isinstance(value, int, bool):
                    integers[name].append((int(value), i))
                else:
                    1 / 0

        self.integers = {}
        for symbol, data in integers.items():
            data.sort()
            ids = []
            indices = [0]
            nmin = data[0][0]
            nmax = data[-1][0]
            assert nmax - nmin < 100, nmax - nmin  # Too wide range!
            m = nmin
            for j, (n, i) in enumerate(data):
                ids.append(i)
                if n > m:
                    indices.append(j)
                    m = n
            indices.append(j + 1)
            self.integers[symbol] = (nmin, nmax, ids, indices)

        self.floats = {}
        for name, data in floats.items():
            assert name not in self.integers
            data.sort()
            ids = [i for value, i in data]
            values = [value for value, i in data]
            self.floats[name] = (ids, values)

    def key(self, name, op, value):
        if name in chemical_symbols:
            n = value
            if name not in self.integers:
                if op == '=' and n != 0:
                    return set()
                if op == '<' and n == 0:
                    return set()
                if op == '>' and n == 0:
                    return set()
                if op == '>=' and n > 0:
                    return set()
                return self.ids
            return self.integer_key(name, op, n)

        if name in self.strings:
            if op == '=':
                if value in self.strings[name]:
                    return self.strings[name][value]
                else:
                    return set()
            if op == '!=':
                result = set()
                for val, ids in self.strings[name].items():
                    if val != value:
                        result.update(ids)
                return result
            1 / 0

        if name in self.floats:
            return self.float_key(name, op, value)

        if name in self.integers:
            return self.integer_key(name, op, value)

        return set()

    def float_key(self, name, op, value):
        values, ids = self.floats[name]
        if op == '!=':
            return self.ids - self.float_key(name, '=', value)
        if op == '<=':
            value = np.nextafter(value, value + 1)
            op = '<'
        if op == '>':
            value = np.nextafter(value, value + 1)
            op = '<='
        j1 = bisect(values, value)
        N = len(values)
        if op == '=':
            result = set()
            while j1 < N and values[j1] == value:
                result.add(ids[j1])
                j1 += 1
            return result
        if op == '<':
            return set(ids[:j1])
        if op == '>=':
            return set(ids[j1:])
        1 / 0

    def integer_key(self, name, op, n):
        nmin, nmax, ids, indices = self.integers[name]
        if op == '=':
            if nmin <= n <= nmax:
                d = n - nmin
                j1, j2 = indices[d:d + 2]
                return set(ids[j1:j2])
            else:
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
            if n > nmax:
                n = nmax
            d = n - nmin
            j = indices[d]
            return set(ids[:j])
        if op == '>':
            if n >= nmax:
                return set()
            if n <= nmin:
                n = nmin
            d = n - nmin
            j = indices[d]
            return set(ids[j:])
        assert False, op

    def formula(self, f):
        ids = None
        for symbol, n in Formula(f).count().items():
            if ids is None:
                ids = self.integer_key(symbol, '>=', n)
            else:
                ids &= self.integer_key(symbol, '>=', n)
        return ids


def bisect(values, value):
    """Find index of value in sorted list of floats.

    >>> bisect([0.0, 1.0, 2.0], 1.5)
    2
    >>> bisect([0.0, 1.0, 2.0], 2.0)
    2
    >>> bisect([0.0, 1.0, 2.0], 2.5)
    3
    """
    if values[0] >= value:
        return 0
    N = len(values)
    if values[-1] < value:
        return N
    j1 = 0
    j2 = N - 1
    while True:
        j = (j1 + j2) // 2
        if values[j] < value:
            j1 = j
        else:
            j2 = j
        if j1 == j2:
            return j


if __name__ == '__main__':
    print(parse1(' '.join(sys.argv[1:])))
