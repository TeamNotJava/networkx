import networkx as nx
import math
from nose.tools import assert_true


def test_graph1():
    G = nx.Graph([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    embedding_data = {0: [1, 2, 3], 1: [2, 0], 2: [3, 0, 1], 3: [2, 0]}
    embedding = nx.PlanarEmbedding()
    embedding.set_data(embedding_data)
    pos = nx.combinatorial_embedding_to_pos(embedding)
    assert_true(is_planar_drawing_correct(G, pos),
                "Planar drawing is not correct")
    assert_true(planar_drawing_conforms_to_embedding(embedding, pos),
                "Planar drawing does not conform to the embedding")


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
    for (a,b) in G.edges():
        for (c,d) in G.edges:
            if(a != c and b != d): #need to have different end points for a chance of conflict
                x1,y1 = pos[a]
                x2,y2 = pos[b]
                x3,y3 = pos[c]
                x4,y4 = pos[d]
                determinant = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
                if(determinant != 0): #the lines are not parallel
                    #calculate intersection point, see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
                    px = (x1*y2 - y1*x2)*(x3 - x4) - (x1-x2)*(x3*y4-y3*x4) /determinant
                    py = (x1*y2 - y1*x2)*(y3 - y4) - (y1-y2)*(x3*y4-y3*x4) / determinant
                    #Check if intersection lies between the points
                    dist_a_b = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
                    dist_a_p= math.sqrt((x1 - px)**2 + (y1 - py)**2)
                    dist_b_p= math.sqrt((x2 - px)**2 + (y2 - py)**2)
                    if(dist_a_p + dist_b_p ==dist_a_b ):
                        print("There is an intersection at {},{}".format(px,py))
                        return False
    return True


def planar_drawing_conforms_to_embedding(embedding, pos):
    """Checks if pos conforms to the planar embedding

    Returns true iff the neighbors are actually oriented in the orientation
    specified of the embedding
    """
    for _,nbrs in embedding:
        if len(nbrs)>1:
            for idx,n in nbrs:
                x1,y1 = pos[n]
                x2,y2 = pos[nbrs[(idx+1)%len(nbrs)]]
                if x1==x2 and y1==y2 : # There should be no nodes mapped to identical positions
                    return False
                if x1 > x2:
                    if y1 < y2: # If the point is right of its predecessor it also has to be lower
                        return False
                if x1 < x2:
                    if y1 > y2:
                        return False
    return True