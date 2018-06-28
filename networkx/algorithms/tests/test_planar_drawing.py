import networkx as nx
from nose.tools import assert_true


def test_graph1():
    G = nx.Graph([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    pos = nx.combinatorial_embedding_to_pos(
        {0: [1, 2, 3], 1: [2, 0], 2: [3, 0, 1], 3: [2, 0]})
    assert_true(is_planar_drawing_correct(G, pos), "Planar drawing is not correct")


def is_planar_drawing_correct(G, pos):
    """Checks if pos represents a planar drawing.

    Check all edges in G for intersections.

    Parameters
    ----------
    G : NetworkX graph
    pos : dict
        Maps every node to a tuple (x, y) representing its position

    Returns
    -------
    is_correct : bool
    """
    pass