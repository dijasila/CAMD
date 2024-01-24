from cxdb.panels.convex_hull import plot_2d
from ase.phasediagram import PhaseDiagram


def test_2d():
    pd = PhaseDiagram(
        [('Mo', 0.0),
          ('Mo2S2', -0.18),
          ('MoS2', -0.92),
          ('S', 0.0)])
    fig = plot_2d(pd)
