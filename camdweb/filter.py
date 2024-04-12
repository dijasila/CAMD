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

from collections import defaultdict

import numpy as np
from ase.data import chemical_symbols
from ase.formula import Formula

from camdweb import ColVal
from camdweb.utils import formula_dict_to_string


class Index:
    def __init__(self,
                 rows: list[tuple[dict[str, int], dict[str, ColVal]]],
                 max_int_range: int = 350,
                 strict=True):
        """Indices for speeding up row filtering.

        The *rows* argument is a list of two element tuples:

        1) Atom count as a dict: {'H': 2, 'O': 1}.

        2) Key-value pairs: {'ok': True, 'gap': 1.1, ...}.

        >>> rows = [({'B': 2, 'N': 2}, {'gap': 3.6}),
        ...         ({'C': 2}, {'gap': 0.0})]
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

        for i, (count, keys) in enumerate(rows):
            nunits = np.gcd.reduce(list(count.values()))
            reduced = formula_dict_to_string(
                {symbol: n // nunits for symbol, n in count.items()})
            self.reduced[reduced].add(i)
            self.natoms[i] = sum(count.values())
            self.ids.add(i)
            for symbol, n in count.items():
                integers[symbol].append((n, i))
            for name, value in keys.items():
                if isinstance(value, str):
                    dct = self.strings[name]
                    if value not in dct:
                        dct[value] = {i}
                    else:
                        dct[value].add(i)
                elif isinstance(value, float):
                    floats[name].append((value, i))
                else:
                    integers[name].append((int(value), i))

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

        sf = self.strings.keys() & self.floats.keys()
        si = self.strings.keys() & self.integers.keys()
        fi = self.floats.keys() & self.integers.keys()
        if sf or si or fi:
            print(sf, si, fi)
            if strict:
                raise ValueError
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
