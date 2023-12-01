import functools
import re
import sys
from typing import Any, Callable

from ase.formula import Formula


@functools.lru_cache
def parse(q: str) -> Callable[[dict[str, int], dict[str, Any]], bool]:
    """Convert query string to Python function.

    >>> f = parse('H2,xc=PBE')
    >>> f({'H': 2}, {'xc': 'PBE'})
    True
    """
    q = parse1(q)
    return eval(f'lambda n, k: {q}')


def parse1(q: str) -> str:
    """Quick'n'dirty hacky parsing of query string to Python expression.

    In the Python expression, *n* is a dict mapping chemical symbols their
    numbers and *k* is a key-value dict:

    >>> parse1('H2')
    "((n.get('H', 0) >= 2))"
    >>> parse1('H2,xc=PBE')
    "((n.get('H', 0) >= 2)) and ('xc' in k and k['xc'] == 'PBE')"
    >>> parse1('H2,Ag=0')
    "((n.get('H', 0) >= 2)) and (n.get('Ag', 0) == 0)"
    >>> parse1('H>7,(Ag=1|Cu=1)')
    "(n.get('H', 0) > 7) and ((n.get('Ag', 0) == 1) or (n.get('Cu', 0) == 1))"
    """
    q = q.replace(' ', '')
    n1 = len(q)
    p1 = q.count('(')
    p2 = q.count(')')
    if p1 != p2:
        raise SyntaxError('Mismatched parenthesis')

    # number of parsed characters (should match n1 when we are done):
    n = p1 + p2

    q = q.replace(',', ' and ')
    n2 = len(q)
    n += (n2 - n1) // 4

    q = q.replace('|', ' or ')
    n3 = len(q)
    n += (n3 - n2) // 3

    h: list[str] = []  # value strings
    for x in ['=', '<', '>', '<=', '>=', '!=']:
        y = '==' if x == '=' else x

        # Find Ag=7, Cu<=7, ...
        matches = reversed(list(re.finditer(fr'([A-Z][a-z]?{x}[0-9]+)', q)))
        for m in matches:
            i, j = m.span()
            k, v = m[1].split(x)
            q = q[:i] + f'(n.get({k!r}, 0) {y} {v})' + q[j:]
            n += j - i

        # Find key=value, key>value, ...
        matches = reversed(
            list(re.finditer(fr'([a-z][a-z0-9_]+{x}[-A-Za-z0-9.+_]+)', q)))
        for m in matches:
            i, j = m.span()
            k, v = m[1].split(x)
            v = repr(str2obj(v))
            q = q[:i] + f'({k!r} in k and k[{k!r}] {y} #{len(h)})' + q[j:]
            h.append(v)
            n += j - i

    # Find Ag7O14, ...
    matches = reversed(list(re.finditer(r'([A-Z][a-z]?[0-9]*)+', q)))
    for m in matches:
        i, j = m.span()
        if j < len(q) and q[j] == "'":
            continue  # we already dealt with this one
        f = Formula(m[0])
        qf = ' and '.join(f'(n.get({s!r}, 0) >= {n})'
                          for s, n in f.count().items())
        q = q[:i] + f'({qf})' + q[j:]
        n += j - i

    # Insert value strings:
    for i, v in reversed(list(enumerate(h))):
        q = q.replace(f'#{i}', v)

    if n1 != n:
        raise SyntaxError('Bad query string')

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


if __name__ == '__main__':
    print(parse1(' '.join(sys.argv[1:])))
