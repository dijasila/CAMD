import pytest
from ase.phasediagram import PhaseDiagram

from camdweb.panels.convex_hull import (group_references,
                                        make_figure_and_tables, plot_2d,
                                        plot_3d, calculate_ehull_energies)


def test_2d():
    pd = PhaseDiagram(
        [('Mo', 0.0),
         ('Mo2S2', -0.18),
         ('MoS2', -0.92),
         ('S', 0.0)])
    plot_2d(pd)


def test_3d():
    pd = PhaseDiagram(
        [('A', 0.0),
         ('ABC', -0.18),
         ('B', 0.0),
         ('C', 0.0)])
    plot_3d(pd)


@pytest.mark.parametrize('n', [1, 2, 3, 4])
def test_1234(n):
    refs = {s: ({s: 1}, 0.0, 'OQMD') for s in 'ABCD'[:n]}
    refs['x'] = ({X: 1 for X in 'ABCD'[:n]}, -0.5 * n, 'C2DB')
    chull, oqmd, c2db = make_figure_and_tables(refs)
    assert '-0.50' in c2db


def test_group_err():
    with pytest.raises(ValueError):
        group_references({'1': ('B', 'A')}, [])


def test_group():
    refs = {'1': ('A',),
            '2': ('B',),
            '3': ('A', 'B'),
            'u1': ('A', 'B'),
            'u2': ('A', 'B')}
    g = group_references(refs, ['u1', 'u2'])
    assert g == {('A', 'B'): ['1', '2', '3', 'u1', 'u2']}


def test_ehull_1d():
    eh = calculate_ehull_energies(
        {'i1': ({'A': 1}, 0.0),
         'i2': ({'A': 2}, 1.0)},
        {'i2'})
    assert eh == {'i2': 1.0}
