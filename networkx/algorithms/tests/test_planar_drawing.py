import networkx as nx
import math
from nose.tools import assert_true


def test_graph1():
    embedding_data = {0: [1, 2, 3], 1: [2, 0], 2: [3, 0, 1], 3: [2, 0]}
    check_embedding_data(embedding_data)


def test_circle_graph():
    G = nx.cycle_graph(20)
    check_embedding_data(G)


def test_grid_graph():
    G = nx.grid_graph(dim=[5, 5])
    check_embedding_data(G)


def test_small_graph():
    embedding_data = {0: [1], 1: [0]}
    check_embedding_data(embedding_data)


def check_embedding_data(embedding_data):
    """ Checks that the planar embedding of the input is correct

    Parameters
    ----------
    embedding_data : NetworkX graph or dict of sorted adjacency lists

    """
    if isinstance(embedding_data, nx.Graph):
        is_planar, embedding = nx.check_planarity(embedding_data)
        assert_true(is_planar, "Graph is not planar")
    else:
        embedding = nx.PlanarEmbedding()
        embedding.set_data(embedding_data)
    pos = nx.combinatorial_embedding_to_pos(embedding)
    assert_true(planar_drawing_conforms_to_embedding(embedding, pos),
                "Planar drawing does not conform to the embedding")
    assert_true(is_planar_drawing_correct(embedding, pos),
                "Intersection in planar drawing")


def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    # Check if float numbers are basically equal, for python >=3.5 there is
    # function for that in the standard library
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def point_inbetween(a, b, p):
    # checks if p is on the line between a and b
    x1, y1 = a
    x2, y2 = b
    px, py = p
    dist_1_2 = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    dist_1_p = math.sqrt((x1 - px)**2 + (y1 - py)**2)
    dist_2_p = math.sqrt((x2 - px)**2 + (y2 - py)**2)
    return is_close(dist_1_p+dist_2_p, dist_1_2)


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
    for (a, b) in G.edges():
        for (c, d) in G.edges():
            if(a != c) and (b != d) and (b != c) and (a != d):   # need to have different end points for conflict
                x1, y1 = pos[a]
                x2, y2 = pos[b]
                x3, y3 = pos[c]
                x4, y4 = pos[d]
                determinant = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
                if determinant != 0:  # the lines are not parallel
                    # calculate intersection point, see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
                    px = (x1*y2 - y1*x2)*(x3 - x4) - (x1-x2)*(x3*y4-y3*x4) / float(determinant)
                    py = (x1*y2 - y1*x2)*(y3 - y4) - (y1-y2)*(x3*y4-y3*x4) / float(determinant)

                    # Check if intersection lies between the points

                    if point_inbetween(pos[a], pos[b], (px, py)) and point_inbetween(pos[c], pos[d], (px, py)):
                        print("There is an intersection at {},{}".format(px, py))
                        return False

                #  Check overlap
                if point_inbetween(pos[a], pos[b], pos[c]) or point_inbetween(pos[a], pos[b], pos[d]):
                    print("A node lies directly on a edge connecting two other nodes")
                    return False
                if point_inbetween(pos[c], pos[d], pos[a]) or point_inbetween(pos[c], pos[d], pos[b]):
                    print("A node lies directly on a edge connecting two other nodes")
                    return False

    return True


class Vector(object):
    """Compare vectors by their angle without loss of precision

    All vectors in direction [0, 1] are the smallest.
    The vectors grow in clockwise direction.
    """
    __slots__ = ['x', 'y', 'node', 'quadrant']

    def __init__(self, x, y, node):
        self.x = x
        self.y = y
        self.node = node
        if self.x >= 0 and self.y > 0:
            self.quadrant = 1
        elif self.x > 0 and self.y <= 0:
            self.quadrant = 2
        elif self.x <= 0 and self.y < 0:
            self.quadrant = 3
        else:
            self.quadrant = 4

    def __eq__(self, other):
        return (self.quadrant == other.quadrant and
                self.x * other.y == self.y * other.x)

    def __lt__(self, other):
        if self.quadrant < other.quadrant:
            return True
        elif self.quadrant > other.quadrant:
            return False
        else:
            return self.x * other.y < self.y * other.x

    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return not other < self

    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return not self < other


def planar_drawing_conforms_to_embedding(embedding, pos):
    """Checks if pos conforms to the planar embedding

    Returns true iff the neighbors are actually oriented in the orientation
    specified of the embedding
    """
    for v in embedding:
        nbr_vectors = []
        v_pos = pos[v]
        for nbr in embedding[v]:
            new_vector = Vector(pos[nbr][0] - v_pos[0], pos[nbr][1] - v_pos[1],
                                nbr)
            nbr_vectors.append(new_vector)
        # Sort neighbors according to their phi angle
        nbr_vectors.sort()
        for idx, nbr_vector in enumerate(nbr_vectors):
            cw_vector = nbr_vectors[(idx + 1) % len(nbr_vectors)]
            ccw_vector = nbr_vectors[idx - 1]
            if (embedding[v][nbr_vector.node]['cw'] != cw_vector.node or
                    embedding[v][nbr_vector.node]['ccw'] != ccw_vector.node):
                return False
            if cw_vector.node != nbr_vector.node and cw_vector == nbr_vector:
                # Lines overlap
                return False
            if ccw_vector.node != nbr_vector.node and ccw_vector == nbr_vector:
                # Lines overlap
                return False
    return True
