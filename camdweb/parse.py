from __future__ import annotations
import sys
import functools
import re
from typing import Callable

from ase.data import chemical_symbols
from camdweb.filter import Index


@functools.lru_cache
def parse(q: str) -> Callable[[Index], set[int]]:
    """Convert filter string to Python function.

    >>> f = parse('H2,xc=PBE')
    >>> i = Index([({'H': 2}, {'xc': 'PBE'})])
    Rows: 1 | Strings: 1 | Integers: 1 | Floats: 0 | Int-floats: 0
    >>> f(i)
    {0}
    >>> f = parse('gap > 5.0')
    >>> i = Index([({'H': 2}, {'gap': 10.0}),
    ...            ({'Si': 2}, {'gap': 1.1})])
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
            v = repr(str2obj(v))
            q = q[:i] + f'(i.key({k!r}, {x!r}, #{len(h)}))' + q[j:]
            h.append(v)
            n += j - i

    # Find Ag7O14, ...
    matches = reversed(list(re.finditer(r'([A-Z][a-z]?[0-9]*)+', q)))
    for m in matches:
        i, j = m.span()
        if j < len(q) and q[j] == "'" or i > 0 and q[i - 1] == "'":
            continue  # we already dealt with this one
        q = q[:i] + f'(i.formula({m[0]!r}))' + q[j:]
        n += j - i

    # Insert value strings:
    for i, v in reversed(list(enumerate(h))):
        q = q.replace(f'#{i}', v)

    if n1 != n:
        raise SyntaxError(f'Bad filter string: {q0!r}')

    print(q)
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


if __name__ == '__main__':
    print(parse1(' '.join(sys.argv[1:])))
