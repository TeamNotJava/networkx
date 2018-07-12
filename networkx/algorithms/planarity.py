from collections import defaultdict
import networkx as nx

__all__ = ["check_planarity"]


def check_planarity(G, counterexample=False):
    """Checks if a graph is planar and returns a counterexample or an embedding

     A graph is planar iff it can be drawn in a plane without
     any edge intersections.

    Parameters
    ----------
    G : NetworkX graph
    counterexample : bool
        A Kuratowski subgraph (to proof non planarity) is only returned if set
        to true

    Returns
    -------
    is_planar : bool
        Is true if the graph is planar

    certificate :
        If the graph is planar this is a planar embedding (dict).
        If the graph is not planar and counterexample is true,
        this is a Kuratowski subgraph.

    Notes
    -----
    A (combinatorial) embedding consists of cyclic orderings of the incident
    edges at each vertex, given such an embedding there are multiple approaches
    discussed in literature to drawing the graph (subject to various
    constraints, e.g. integer coordinates), see e.g. [2].

    The planarity check algorithm and extraction of the combinatorial embedding
    is based on the Left-Right Planarity Test [1].

    A counterexample is only generated if the corresponding parameter is set,
    because the complexity of the counterexample generation is higher.

    References
    ----------
    .. [1] Ulrik Brandes:
        The Left-Right Planarity Test
        2009
        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.217.9208
    .. [2] Takao Nishizeki, Md Saidur Rahman:
        Planar graph drawing
        Lecture Notes Series on Computing: Volume 12
        2004
    """

    planarity_state = LRPlanarity(G)
    embedding = planarity_state.lr_planarity()
    if embedding is None:
        # graph is not planar
        if counterexample:
            return False, get_counterexample(G)
        else:
            return False, None
    else:
        # graph is planar
        return True, embedding.get_data()


def get_counterexample(G):
    """Obtains a Kuratowski subgraph

    Raises nx.NetworkXException if G is planar.

    The function removes edges such that the graph is still not planar.
    At some point the removal of any edge would make the graph planar.
    This subgraph must be a Kuratowski subgraph.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    subgraph : NetworkX graph
        A Kuratowski subgraph that proves that G is not planar.

    """
    # copy graph
    G = nx.Graph(G)

    if check_planarity(G)[0]:
        raise nx.NetworkXException("G is planar - no counter example.")

    # find Kuratowski subgraph
    subgraph = nx.Graph()
    for u in G:
        nbrs = list(G[u])
        for v in nbrs:
            G.remove_edge(u, v)
            if check_planarity(G)[0]:
                G.add_edge(u, v)
                subgraph.add_edge(u, v)

    return subgraph


class Interval(object):
    """Represents a set of return edges

    All return edges in an interval induce a same constraint on the contained
    edges, which means that all edges must either have a left orientation or
    all edges must have a right orientation.
    """
    def __init__(self, low=None, high=None):
        self.low = low
        self.high = high

    def empty(self):
        """Check if the interval is empty"""
        return self.low is None and self.high is None

    def copy(self):
        """Return a copy of this interval"""
        return Interval(self.low, self.high)

    def conflicting(self, b, planarity_state):
        """Return True if interval I conflicts with edge b"""
        return (not self.empty() and
                planarity_state.lowpt[self.high] > planarity_state.lowpt[b])


class ConflictPair(object):
    """Represents a different constraint between two intervals

    The edges in the left interval must have a different orientation than
    the one in the right interval.
    """
    def __init__(self, left=Interval(), right=Interval()):
        self.left = left
        self.right = right

    def swap(self):
        """Swap left and right intervals"""
        temp = self.left
        self.left = self.right
        self.right = temp

    def lowest(self, planarity_state):
        """Return the lowest lowpoint of a conflict pair"""
        if self.left.empty():
            return planarity_state.lowpt[self.right.low]
        if self.right.empty():
            return planarity_state.lowpt[self.left.low]
        return min(planarity_state.lowpt[self.left.low],
                   planarity_state.lowpt[self.right.low])


def top_of_stack(l):
    """Returns the element on top of the stack."""
    if not l:
        return None
    return l[-1]


class LRPlanarity(object):
    """A class to maintain the state during planarity check"""
    def __init__(self, G):
        # copy G without adding self-loops
        self.G = nx.Graph()
        self.G.add_nodes_from(G.nodes)
        for e in G.edges:
            if e[0] != e[1]:
                self.G.add_edge(e[0], e[1])

        self.roots = []

        # distance from tree root
        self.height = defaultdict(lambda: None)

        self.lowpt = {}  # height of lowest return point of an edge
        self.lowpt2 = {}  # height of second lowest return point
        self.nesting_depth = {}  # for nesting order

        # None -> missing edge
        self.parent_edge = defaultdict(lambda: None)

        # oriented DFS graph
        self.DG = nx.DiGraph()
        self.DG.add_nodes_from(G.nodes)

        self.ordered_adjs = {}

        self.ref = defaultdict(lambda: None)
        self.side = defaultdict(lambda: 1)

        # stack of conflict pairs
        self.S = []
        self.stack_bottom = {}
        self.lowpt_edge = {}

        self.left_ref = {}
        self.right_ref = {}

        self.embedding = PlanarEmbedding()

    def lr_planarity(self):
        """Execute the LR planarity test

        Returns
        -------
        embedding: dict
            If the graph is planar an embedding is returned. Otherwise None.
        """
        if self.G.order() > 2 and self.G.size() > 3 * self.G.order() - 6:
            # graph is not planar
            return None

        # orientation of the graph by depth first search traversal
        for v in self.G:
            if self.height[v] is None:
                self.height[v] = 0
                self.roots.append(v)
                self.dfs_orientation(v)

        self.G = None  # just unsetting this for correctness purposes

        # testing
        for v in self.DG:  # sort the adjacency lists by nesting depth
            # note: this sorting leads to non linear time
            self.ordered_adjs[v] = sorted(
                self.DG[v], key=lambda w: self.nesting_depth[(v, w)])
        for v in self.roots:
            if not self.dfs_testing(v):
                return None

        for e in self.DG.edges:
            self.nesting_depth[e] = self.sign(e) * self.nesting_depth[e]

        self.embedding.add_nodes_from(self.DG.nodes)
        for v in self.DG:
            # sort the adjacency lists again
            self.ordered_adjs[v] = sorted(
                self.DG[v], key=lambda w: self.nesting_depth[(v, w)])
            # initialize the embedding
            previous_node = None
            for w in self.ordered_adjs[v]:
                self.embedding.add_half_edge_cw(v, w, previous_node)
                previous_node = w

        # compute the complete embedding
        for v in self.roots:
            self.dfs_embedding(v)

        return self.embedding

    # orient the graph by DFS-traversal, compute lowpoints and nesting order
    def dfs_orientation(self, v):
        e = self.parent_edge[v]
        for w in self.G[v]:
            if (v, w) in self.DG.edges or (w, v) in self.DG.edges:
                continue  # the edge was already oriented
            vw = (v, w)
            self.DG.add_edge(v, w)  # orient the edge

            self.lowpt[vw] = self.height[v]
            self.lowpt2[vw] = self.height[v]
            if self.height[w] is None:  # (v, w) is a tree edge
                self.parent_edge[w] = vw
                self.height[w] = self.height[v] + 1
                self.dfs_orientation(w)
            else:  # (v, w) is a back edge
                self.lowpt[vw] = self.height[w]

            # determine nesting graph
            self.nesting_depth[vw] = 2 * self.lowpt[vw]
            if self.lowpt2[vw] < self.height[v]:  # chordal
                self.nesting_depth[vw] += 1

            # update lowpoints of parent edge e
            if e is not None:
                if self.lowpt[vw] < self.lowpt[e]:
                    self.lowpt2[e] = min(self.lowpt[e], self.lowpt2[vw])
                    self.lowpt[e] = self.lowpt[vw]
                elif self.lowpt[vw] > self.lowpt[e]:
                    self.lowpt2[e] = min(self.lowpt2[e], self.lowpt[vw])
                else:
                    self.lowpt2[e] = min(self.lowpt2[e], self.lowpt2[vw])

    def dfs_testing(self, v):
        """Test for LR partition"""
        e = self.parent_edge[v]
        for w in self.ordered_adjs[v]:
            ei = (v, w)
            self.stack_bottom[ei] = top_of_stack(self.S)
            if ei == self.parent_edge[w]:  # tree edge
                if not self.dfs_testing(w):
                    return False
            else:  # back edge
                self.lowpt_edge[ei] = ei
                self.S.append(ConflictPair(right=Interval(ei, ei)))

            # integrate new return edges
            if self.lowpt[ei] < self.height[v]:
                if w == self.ordered_adjs[v][0]:  # e_i has return edge
                    self.lowpt_edge[e] = self.lowpt_edge[ei]
                else:  # add constraints of e_i
                    if not self.add_constraints(ei, e):
                        # graph is not planar
                        return False

        # remove back edges returning to parent
        if e is not None:  # v isn't root
            self.remove_back_edges(e)
        return True

    def add_constraints(self, ei, e):
        P = ConflictPair()
        # merge return edges of e_i into P.right
        while True:
            Q = self.S.pop()
            if not Q.left.empty():
                Q.swap()
            if not Q.left.empty():  # not planar
                return False
            if self.lowpt[Q.right.low] > self.lowpt[e]:
                # merge intervals
                if P.right.empty():  # topmost interval
                    P.right = Q.right.copy()
                else:
                    self.ref[P.right.low] = Q.right.high
                P.right.low = Q.right.low
            else:  # align
                self.ref[Q.right.low] = self.lowpt_edge[e]
            if top_of_stack(self.S) == self.stack_bottom[ei]:
                break
        # merge conflicting return edges of e_1,...,e_i-1 into P.L
        while (top_of_stack(self.S).left.conflicting(ei, self) or
               top_of_stack(self.S).right.conflicting(ei, self)):
            Q = self.S.pop()
            if Q.right.conflicting(ei, self):
                Q.swap()
            if Q.right.conflicting(ei, self):  # not planar
                return False
            # merge interval below lowpt(e_i) into P.R
            self.ref[P.right.low] = Q.right.high
            if Q.right.low is not None:
                P.right.low = Q.right.low

            if P.left.empty():  # topmost interval
                P.left = Q.left.copy()
            else:
                self.ref[P.left.low] = Q.left.high
            P.left.low = Q.left.low

        if not (P.left.empty() and P.right.empty()):
            self.S.append(P)
        return True

    def remove_back_edges(self, e):
        u = e[0]
        # trim back edges ending at parent u
        # drop entire conflict pairs
        while self.S and top_of_stack(self.S).lowest(self) == self.height[u]:
            P = self.S.pop()
            if P.left.low is not None:
                self.side[P.left.low] = -1

        if self.S:  # one more conflict pair to consider
            P = self.S.pop()
            # trim left interval
            while P.left.high is not None and P.left.high[1] == u:
                P.left.high = self.ref[P.left.high]
            if P.left.high is None and P.left.low is not None:
                # just emptied
                self.ref[P.left.low] = P.right.low
                self.side[P.left.low] = -1
                P.left.low = None
            # trim right interval
            while P.right.high is not None and P.right.high[1] == u:
                P.right.high = self.ref[P.right.high]
            if P.right.high is None and P.right.low is not None:
                # just emptied
                self.ref[P.right.low] = P.left.low
                self.side[P.right.low] = -1
                P.right.low = None
            self.S.append(P)

        # side of e is side of a highest return edge
        if self.lowpt[e] < self.height[u]:  # e has return edge
            hl = top_of_stack(self.S).left.high
            hr = top_of_stack(self.S).right.high

            if hl is not None and (
                    hr is None or self.lowpt[hl] > self.lowpt[hr]):
                self.ref[e] = hl
            else:
                self.ref[e] = hr

    # complete the embedding
    def dfs_embedding(self, v):
        for w in self.ordered_adjs[v]:
            ei = (v, w)
            if ei == self.parent_edge[w]:  # tree edge
                self.embedding.add_half_edge_first(w, v)
                self.left_ref[v] = w
                self.right_ref[v] = w
                self.dfs_embedding(w)
            else:  # back edge
                if self.side[ei] == 1:
                    # place v directly after right_ref[w] in embed. list of w
                    self.embedding.add_half_edge_cw(w, v, self.right_ref[w])
                else:
                    # place v directly before left_ref[w] in embed. list of w
                    self.embedding.add_half_edge_ccw(w, v, self.left_ref[w])
                    self.left_ref[w] = v

    # function to resolve the relative side of an edge to the absolute side
    def sign(self, e):
        if self.ref[e] is not None:
            self.side[e] = self.side[e] * self.sign(self.ref[e])
            self.ref[e] = None
        return self.side[e]


class PlanarEmbedding:
    """ Represents a planar embedding.

    This class maintains an order on the outgoing edges for each node.
    In order to represent a valid planar embedding of a Graph each edge of the
    graph must be represented by two half edges in either direction.

    There are three different ways to add edges to the planar embedding:
    - add_half_edge_{ccw,cw} : Inserts one half edge at a specific position.
        This method takes constant amount of time.
        The resulting embedding is not guaranteed to be valid.

    - add_edge_{ccw,cw} : Inserts two half edges at a specific position.
        This method cannot be used if the nodes are in different components
        The method can takes more time for larger graphs.
        If no exception is thrown the embedding stays valid afterwards.

    - add_edge : Adds two half edges between some random edges.
        This method takes a constant amount of time.
        The resulting embedding is only guaranteed to be valid, if the
        connected nodes were contained in different graph components.
    """

    def __init__(self):
        # Maps nodes to a dict mapping neighbor nodes to the ccw neighbor node
        self.ccw_nbr = {}
        # Maps nodes to a dict mapping neighbor nodes to the cw neighbor node
        self.cw_nbr = {}
        # Maps a node to the first neighbor
        self.first_nbr = {}

    def get_data(self):
        """Converts this object into a dict of list of nodes structure"""
        embedding = dict()
        for v in self.ccw_nbr:
            embedding[v] = list(self.get_neighbors(v))
        return embedding

    def get_neighbors(self, v):
        """Yields the neighbors of v in clockwise order"""
        if v not in self.first_nbr:
            # v has no neighbors
            return
        start_node = self.first_nbr[v]
        nbr_order = self.ccw_nbr[v]
        yield start_node
        current_node = nbr_order[start_node]
        while start_node != current_node:
            yield current_node
            current_node = nbr_order[current_node]

    def check_structure(self):
        """Returns true if every half edge has its opposite half edge"""
        pass

    def check_intersection(self):
        """Returns true if there are no intersections of edges"""
        pass

    def add_node(self, v):
        """Adds a node to the embedding if it does not exist"""
        if v not in self.ccw_nbr:
            self.ccw_nbr[v] = {}
            self.cw_nbr[v] = {}

    def add_nodes_from(self, nodes):
        for v in nodes:
            self.add_node(v)

    def add_half_edge_ccw(self, start_node, end_node, reference_neighbor, update_fist_nbr=True):
        """Adds a half edge from start_node to end_node.

        The half edge is added counter clockwise next to the existing half edge
        (start_node, reference_neighbor).

        Calling this method can break the embedding.

        Raises an exception if the reference half edge does not exist.

        If there are no hash table collisions the complexity is constant.
        """
        ccw_order = self.ccw_nbr[start_node]
        cw_order = self.cw_nbr[start_node]
        if reference_neighbor is None and len(ccw_order) == 0:
            # The start node has no neighbors
            ccw_order[end_node] = end_node
            cw_order[end_node] = end_node
            self.first_nbr[start_node] = end_node
            return
        if reference_neighbor not in ccw_order:
            raise nx.NetworkXException(
                "Cannot add edge. Reference neighbor does not exist")
        # Get half edge at the other side
        ccw_reference = ccw_order[reference_neighbor]
        # Alter half edge data structures
        ccw_order[reference_neighbor] = end_node
        ccw_order[end_node] = ccw_reference
        cw_order[ccw_reference] = end_node
        cw_order[end_node] = reference_neighbor

        if update_fist_nbr and reference_neighbor == self.first_nbr[start_node]:
            # Update first neighbor
            self.first_nbr[start_node] = end_node

    def add_half_edge_first(self, start_node, end_node):
        reference = self.first_nbr.get(start_node, None)
        self.add_half_edge_ccw(start_node, end_node, reference)

    def add_half_edge_cw(self, start_node, end_node, reference_neighbor):
        """Adds a half edge from start_node to end_node.

        The half edge is added clockwise next to the existing half edge
        (start_node, reference_neighbor).

        Calling this method can break the embedding.

        Raises an exception if the reference half edge does not exist, or if
        adding the specified edge would break the planar embedding.

        If there are no hash table collisions the complexity is constant.
        """

        if reference_neighbor is None and len(self.cw_nbr[start_node]) == 0:
            # The start node has no neighbors
            self.ccw_nbr[start_node][end_node] = end_node
            self.cw_nbr[start_node][end_node] = end_node
            self.first_nbr[start_node] = end_node
        else:
            cw_reference = self.cw_nbr[start_node][reference_neighbor]
            self.add_half_edge_ccw(start_node, end_node, cw_reference, False)

    def add_edge_ccw(self, start_node, end_node, reference_neighbor):
        """Adds the half edges from start_node to end_node and the reverse.

        The half edge (start_node, end_node) is added counter clockwise next to
        the
        existing half edge (start_node, reference_neighbor). The half edge
        in the other direction is automatically determined by traversing the
        face that the added edge splits.

        Calling this method on a valid planar embedding object guarantees that
        the resulting object is still a valid planar embedding.

        Raises an exception if the reference half edge does not exist, if
        adding the specified edge would break the planar embedding or if
        start_node and end_node are in different components.

        The complexity is linear in the size of the graph.
        """
        if reference_neighbor not in self.ccw_nbr[start_node]:
            nx.NetworkXException(
                "Cannot add edge. The start node has degree 0.")
        ccw_reference = self.ccw_nbr[start_node][reference_neighbor]
        self.add_edge_ccw(start_node, end_node, ccw_reference)

    def add_edge_cw(self, start_node, end_node, reference_neighbor):
        """Adds the half edges from start_node to end_node and the reverse.

        The half edge (start_node, end_node) is added clockwise next to
        the existing half edge (start_node, reference_neighbor). The half edge
        in the other direction is automatically determined by traversing the
        face that the added edge splits.

        Calling this method on a valid planar embedding object guarantees that
        the resulting object is still a valid planar embedding.

        Raises an exception if the reference half edge does not exist, if
        adding the specified edge would break the planar embedding or if
        start_node and end_node are in different components.

        The complexity is linear in the size of the graph.
        """
        prev_face_node = start_node
        for face_node in self.traverse_face(start_node, reference_neighbor):
            if face_node == end_node:
                self.add_half_edge_cw(start_node, end_node, reference_neighbor)
                self.add_half_edge_ccw(end_node, start_node, prev_face_node)
                return
            prev_face_node = face_node
        exception_msg = "Cannot add edge. End node is either in a different " \
                        "component, or the edge would would violate planarity."
        raise nx.NetworkXException(exception_msg)

    def add_edge(self, v, w):
        """Adds half edges for (v, w) and (w, v) at some position.

        This method should only be called if v and w are in different
        components, or it might break the embedding.

        The complexity is constant.
        """
        if self.has_edge(v, w):
            # Edge already present
            raise nx.NetworkXException("Edge is already present")

        if len(self.cw_nbr[v]) == 0:
            self.cw_nbr[v][w] = w
            self.ccw_nbr[v][w] = w
            self.first_nbr[v] = w
        else:
            reference = next(iter(self.cw_nbr[v]))
            self.add_half_edge_ccw(v, w, reference)

        if len(self.cw_nbr[w]) == 0:
            self.cw_nbr[w][v] = v
            self.ccw_nbr[w][v] = v
            self.first_nbr[w] = v
        else:
            reference = next(iter(self.cw_nbr[w]))
            self.add_half_edge_ccw(w, v, reference)

    def has_edge(self, v, w):
        """Returns true if both half edges (v, w) and (w, v) are present"""
        return w in self.ccw_nbr[v]

    def traverse_face(self, v, w, mark_half_edges=None):
        """Yields nodes on the face of the half edge (v, w)

        The face lies to the right of the half edge (when v is at the bottom
        and w at the top).

        Optionally it is possible to pass a set in which all encountered half
        edges are added.
        """
        yield v
        if mark_half_edges is not None:
            mark_half_edges.add((v, w))
        previous_node = v
        current_node = w
        while current_node != v:
            yield current_node
            previous_node = current_node
            current_node = self.ccw_nbr[current_node][previous_node]
            if mark_half_edges is not None:
                mark_half_edges.add((previous_node, current_node))

    def remove_half_edge(self, v, w):
        """ Removes the half edge (v, w)
        """
        if len(self.ccw_nbr[v]) == 1:
            self.ccw_nbr[v] = {}
            del self.first_nbr[v]
        else:
            ccw_neighbor = self.ccw_nbr[v][w]
            cw_neighbor = self.cw_nbr[v][w]
            self.cw_nbr[v][ccw_neighbor] = cw_neighbor
            self.ccw_nbr[v][cw_neighbor] = ccw_neighbor
            self.first_nbr[v] = cw_neighbor

    def remove_edge(self, v, w):
        """Removes both half edges between v and w.

        Raises an exception if not BOTH half edges are present.
        The complexity is constant."""
        if not self.has_edge(v, w):
            raise nx.NetworkXException("Cannot remove edge. Edge not present.")
        self.remove_half_edge(v, w)
        self.remove_half_edge(w, v)
