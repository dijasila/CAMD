from __future__ import annotations

import functools
import operator
from typing import Any, Callable

import numpy as np
from ase.data import atomic_numbers
from ase.formula import Formula
from lark.lexer import Token  # type: ignore
from lark.tree import Tree  # type: ignore

from optimade.filterparser import LarkParser  # type: ignore
from camdweb.filter import Index

OPS = {
    '=': operator.eq,
    '!=': operator.ne,
    '>': operator.gt,
    '<': operator.lt,
    '>=': operator.ge,
    '<=': operator.le}

REVERSE_OPS = {
    '>': '<',
    '<': '>',
    '>=': '<=',
    '<=': '>='}

SPECIAL = {
    'nsites',
    'nelements',
    'nperiodic_dimensions',
    'chemical_formula_descriptive',
    'chemical_formula_reduced',
    'chemical_formula_anonymous',
    'chemical_formula_hill'}

LENGTH_ALIASES = {
    'elements': 'nelements',
    'element_ratios': 'nelements',
    'cartesian_site_positions': 'nsites',
    'species_at_sites': 'nsites',
}


def select(node: Any,
           index: Index):
    """Create SELECT SQL-statement."""
    if len(node) == 3:
        n1, n2, n3 = node
        if n2[0] == 'OPERATOR' and n3[0] == 'IDENTIFIER':
            key = n3[1]
            op = n2[1]
            value = n1
            op = REVERSE_OPS.get(op, op)
            return index.key(key, op, value)
        raise ValueError

    n1, n2 = node
    if n1 == 'OR':
        v1, v2 = n2
        return select(v1, index) | select(v2, index)
    if n1 == 'AND':
        v1, v2 = n2
        return select(v1, index) & select(v2, index)
    if n1 == ('NOT', 'NOT'):
        return index.ids - select(n2, index)
    if n1[0] == 'IDENTIFIER':
        key = n1[1]
        *n3, n4 = n2
        name = ' '.join(n[0] for n in n3)
        if name.startswith('HAS'):
            if key not in ['elements', 'species_at_sites']:
                raise ValueError
            if name == 'HAS':
                value = n4
                return index.key(value, '>', 0)
            if name == 'HAS ALL':
                values = n4 if isinstance(n4, list) else [n4]
                return functools.reduce(
                    (operator.and_,
                     index.key(value, '>', 0) for value in values))
            assert name == 'HAS ANY', name
            values = n4 if isinstance(n4, list) else [n4]
            return functools.reduce(
                (operator.or_,
                 index.key(value, '>', 0) for value in values))
        if name.startswith('LENGTH'):
            value = n4
            op = "="
            if name.endswith('OPERATOR'):
                op = n3[1][1]
            if key in LENGTH_ALIASES:
                return self.key(LENGTH_ALIASES[key], op, value)
            raise NotImplementedError(
            f'Length filter not supported on field {key!r}')
        if n3[0][0] == 'OPERATOR':
            op = n3[0][1]
            value = n4
            return self.compare(key, op, value)


def parse_lark_tree(node: Tree | Token) -> Any:
    """Convert Lark tree to simple data structure.

    See examples in ``parser_test.py``.
    """
    if isinstance(node, Token):
        if node.type == 'SIGNED_INT':
            return int(node.value)
        if node.type == 'SIGNED_FLOAT':
            return float(node.value)
        if node.type == 'ESCAPED_STRING':
            return node.value[1:-1]
        return (node.type, node.value)
    children = [parse_lark_tree(child)  # type: ignore
                for child in node.children]
    if len(children) == 1:
        return children[0]
    if node.data == 'expression':
        return ('OR', children)
    if node.data == 'expression_clause':
        return ('AND', children)
    return children


def parse(filter: str, parser: LarkParser) -> Any:
    """Parse OPTIMADE filter string to simple data structure"""
    return parse_lark_tree(parser.parse(filter))


def create_parse_function() -> Callable[[str], Any]:
    """Create a parser function."""
    parser = LarkParser(version=(1, 0, 0))
    return functools.partial(parse, parser=parser)
