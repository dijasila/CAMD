import pytest
from cxdb.panels.convex_hull import plot_2d, plot_3d, make_figure_and_tables
from ase.phasediagram import PhaseDiagram


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
