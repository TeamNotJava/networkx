#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Author:  Aric Hagberg (hagberg@lanl.gov),
#          Pieter Swart (swart@lanl.gov),
#          Dan Schult(dschult@colgate.edu)
"""Base class for undirected maps.

The Map class allows any hashable object as a node
and can associate key/value attribute pairs with each undirected edge.

Self-loops are allowed but multiple edges are not (see MultiGraph).

For directed graphs see DiGraph and MultiDiGraph.
"""
from __future__ import division
import warnings
from copy import deepcopy
from collections import Mapping

import networkx as nx
from networkx.classes import Graph
from networkx.classes.coreviews import AtlasView, AdjacencyView
from networkx.classes.reportviews import NodeView, EdgeView, DegreeView
from networkx.exception import NetworkXError
import networkx.convert as convert
from networkx.utils import pairwise


class Map(OrderedGraph):
    """
    Base class for undirected maps.

    A Map stores nodes and edges with optional data, or attributes and the cw order of the edges.

    Graphs hold undirected edges.  Self loops are allowed but multiple
    (parallel) edges are not.

    Nodes can be arbitrary (hashable) Python objects with optional
    key/value attributes. By convention `None` is not used as a node.

    Edges are represented as links between nodes with optional
    key/value attributes.

    Parameters
    ----------
    incoming_graph_data : input graph (optional, default: None)
        Data to initialize graph. If None (default) an empty
        graph is created.  The data can be any format that is supported
        by the to_networkx_graph() function, currently including edge list,
        dict of dicts, dict of lists, NetworkX graph, NumPy matrix
        or 2d ndarray, SciPy sparse matrix, or PyGraphviz graph.

    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to graph as key=value pairs.

    See Also
    --------
    Graph
    """

    def __init__(self, incoming_graph_data=None, **attr):
        """Initialize a map with edges, name, or graph attributes.

        Parameters
        ----------
        incoming_graph_data : input graph
            Data to initialize graph.  If incoming_graph_data=None (default)
            an empty graph is created.  The data can be an edge list, or any
            NetworkX graph object.  If the corresponding optional Python
            packages are installed the data can also be a NumPy matrix
            or 2d ndarray, a SciPy sparse matrix, or a PyGraphviz graph.

        attr : keyword arguments, optional (default= no attributes)
            Attributes to add to graph as key=value pairs.

        See Also
        --------
        Graph
        """
        self._map = dict()
        Graph.__init__(self, incoming_graph_data, **attr)
        # todo: Handle self._map when incoming_graph_Data != None

    def add_edge(self, u_of_edge, , v_of_edge, prev_edge_of_u=None, prev_edge_of_v=None, **attr):
        """Add an edge between u and v.

        The nodes u and v will be automatically added if they are
        not already in the graph.

        Edge attributes can be specified with keywords or by directly
        accessing the edge's attribute dictionary. See examples below.

        Parameters
        ----------
        u, v : nodes
            Nodes can be, for example, strings or numbers.
            Nodes must be hashable (and not None) Python objects.
        prev_u, prev_v : nodes
            An edge with (u, prev_u) must be present. 
            At node u the new edge (u, v) will be added in cw order to node u.
            Respectively for v
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edges_from : add a collection of edges

        Notes
        -----
        Adding an edge that already exists updates the edge data.

        Many NetworkX algorithms designed for weighted graphs use
        an edge attribute (by default `weight`) to hold a numerical value.

        Examples
        --------
        The following all add the edge e=(1, 2) to graph G:

        >>> G = nx.Map()
        >>> G.add_edge(1, None, 2, None)           # explicit two-node form
        >>> G.add_edge(1, None, 3, None)           # explicit two-node form
        >>> G.add_edge(2, None, 3, None)           # explicit two-node form
        >>> G.add_edge(4, None, 1, 3)           # explicit two-node form
        >>> G.add_edge(4, None, 2, 1)           # explicit two-node form
        >>> G.add_edge(4, None, 3, 2)           # explicit two-node form
        """
        u, v = u_of_edge, v_of_edge
        # add nodes
        if u not in self._node:
            self._adj[u] = self.adjlist_inner_dict_factory()
            self._map[u] = list()
            self._node[u] = {}
        if v not in self._node:
            self._adj[v] = self.adjlist_inner_dict_factory()
            self._map[v] = list()
            self._node[v] = {}
        # add the edge
        datadict = self._adj[u].get(v, self.edge_attr_dict_factory())
        datadict.update(attr)
        # Still filling in the adjecency matrix for speed ups
        self._adj[u][v] = datadict
        self._adj[v][u] = datadict

        
        if len(self._map[u]) is 0:
            self._map[u] = [v]  
        else:
            if prev_edge_of_u is not None:
                if prev_edge_of_u not in self._map[u]:
                    raise nx.exception.NodeNotFound('Given previous node {0} not present at node {1}'.format(prev_edge_of_u, u))
                i = self._map[u].index(prev_edge_of_u)
                self._map[u].insert(i+1, v)
            else:
                self._map[u].append(v)

        if len(self._map[v]) is 0:
            self._map[v] = [u]  
        else:
            if prev_edge_of_v is not None:
                if prev_edge_of_v not in self._map[v]:
                    raise nx.exception.NodeNotFound('Given previous node {0} not present at node {1}'.format(prev_edge_of_v, v))
                i = self._map[v].index(prev_edge_of_v)
                self._map[v].insert(i+1, u)
            else:
                self._map[v].append(u)

    def check(self):
        pass

    def add_edges_from(self, ebunch_to_add, **attr):
        """Add all the edges in ebunch_to_add.

        Parameters
        ----------
        ebunch_to_add : container of edges
            Each edge given in the container will be added to the
            graph. The edges must be given as as 2-tuples (u, v) or
            3-tuples (u, v, d) where d is a dictionary containing edge data.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edge : add a single edge
        add_weighted_edges_from : convenient way to add weighted edges

        Notes
        -----
        Adding the same edge twice has no effect but any edge data
        will be updated when each duplicate edge is added.

        Edge attributes specified in an ebunch take precedence over
        attributes specified via keyword arguments.

        Examples
        --------
        >>> G = nx.Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edges_from([(0, 1), (1, 2)]) # using a list of edge tuples
        >>> e = zip(range(0, 3), range(1, 4))
        >>> G.add_edges_from(e) # Add the path graph 0-1-2-3

        Associate data to edges

        >>> G.add_edges_from([(1, 2), (2, 3)], weight=3)
        >>> G.add_edges_from([(3, 4), (1, 4)], label='WN2898')
        """
        for e in ebunch_to_add:
            ne = len(e)
            if ne == 5:
                u, v, prev_edge_of_u, prev_edge_of_v, dd = e
            elif ne == 4:
                u, v, prev_edge_of_u, prev_edge_of_v = e
                dd = {}  # doesn't need edge_attr_dict_factory
            else:
                raise NetworkXError(
                    "Edge tuple %s must be a 4-tuple or 5-tuple." % (e,))
            if u not in self._node:
                self._adj[u] = self.adjlist_inner_dict_factory()
                self._map[u] = list()
                self._node[u] = {}
            if v not in self._node:
                self._adj[v] = self.adjlist_inner_dict_factory()
                self._map[v] = list()
                self._node[v] = {}
            datadict = self._adj[u].get(v, self.edge_attr_dict_factory())
            datadict.update(attr)
            datadict.update(dd)
            self._adj[u][v] = datadict
            self._adj[v][u] = datadict

            if len(self._map[u]) is 0:
                self._map[u] = [v]  
            else:
                if prev_edge_of_u is not None:
                    if prev_edge_of_u not in self._map[u]:
                        raise nx.exception.NodeNotFound('Given previous node {0} not present at node {1}'.format(prev_edge_of_u, u))
                    i = self._map[u].index(prev_edge_of_u)
                    self._map[u].insert(i+1, v)
                else:
                    self._map[u].append(v)

            if len(self._map[v]) is 0:
                self._map[v] = [u]  
            else:
                if prev_edge_of_v is not None:
                    if prev_edge_of_v not in self._map[v]:
                        raise nx.exception.NodeNotFound('Given previous node {0} not present at node {1}'.format(prev_edge_of_v, v))
                    i = self._map[v].index(prev_edge_of_v)
                    self._map[v].insert(i+1, u)
                else:
                    self._map[v].append(u)

    def remove_edge(self, u, v):
        """Remove the edge between u and v.

        Parameters
        ----------
        u, v : nodes
            Remove the edge between nodes u and v.

        Raises
        ------
        NetworkXError
            If there is not an edge between u and v.

        See Also
        --------
        remove_edges_from : remove a collection of edges

        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, etc
        >>> G.remove_edge(0, 1)
        >>> e = (1, 2)
        >>> G.remove_edge(*e) # unpacks e from an edge tuple
        >>> e = (2, 3, {'weight':7}) # an edge with attribute data
        >>> G.remove_edge(*e[:2]) # select first part of edge tuple
        """
        try:
            del self._adj[u][v]
            if u != v:  # self-loop needs only one entry removed
                del self._adj[v][u]
            self._map[u].remove(v)
            self._map[v].remove(u)
        except KeyError:
            raise NetworkXError("The edge %s-%s is not in the graph" % (u, v))

    def remove_edges_from(self, ebunch):
        """Remove all edges specified in ebunch.

        Parameters
        ----------
        ebunch: list or container of edge tuples
            Each edge given in the list or container will be removed
            from the graph. The edges can be:

                - 2-tuples (u, v) edge between u and v.
                - 3-tuples (u, v, k) where k is ignored.

        See Also
        --------
        remove_edge : remove a single edge

        Notes
        -----
        Will fail silently if an edge in ebunch is not in the graph.

        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> ebunch=[(1, 2), (2, 3)]
        >>> G.remove_edges_from(ebunch)
        """
        adj = self._adj
        for e in ebunch:
            u, v = e[:2]  # ignore edge data if present
            if u in adj and v in adj[u]:
                del adj[u][v]
                if u != v:  # self loop needs only one entry removed
                    del adj[v][u]
                self._map[u].remove(v)
                self._map[v].remove(u)

    def clear(self):
        """Remove all nodes and edges from the graph.

        This also removes the name, and all graph, node, map, and edge attributes.

        See Also
        --------
        Graph.clear

        """
        self._adj.clear()
        self._node.clear()
        self._map.clear()
        self.graph.clear()

    def fresh_copy(self):
        """Return a fresh copy graph with the same data structure.

        A fresh copy has no nodes, edges or graph attributes. It is
        the same data structure as the current graph. This method is
        typically used to create an empty version of the graph.

        Notes
        -----
        If you subclass the base class you should overwrite this method
        to return your class of graph.
        """
        return Map()

    def copy(self, as_view=False):
        """Return a copy of the graph.

        The copy method by default returns a shallow copy of the graph
        and attributes. That is, if an attribute is a container, that
        container is shared by the original an the copy.
        Use Python's `copy.deepcopy` for new containers.

        If `as_view` is True then a view is returned instead of a copy.

        Notes
        -----
        All copies reproduce the graph structure, but data attributes
        may be handled in different ways. There are four types of copies
        of a graph that people might want.

        Deepcopy -- The default behavior is a "deepcopy" where the graph
        structure as well as all data attributes and any objects they might
        contain are copied. The entire graph object is new so that changes
        in the copy do not affect the original object. (see Python's
        copy.deepcopy)

        Data Reference (Shallow) -- For a shallow copy the graph structure
        is copied but the edge, node and graph attribute dicts are
        references to those in the original graph. This saves
        time and memory but could cause confusion if you change an attribute
        in one graph and it changes the attribute in the other.
        NetworkX does not provide this level of shallow copy.

        Independent Shallow -- This copy creates new independent attribute
        dicts and then does a shallow copy of the attributes. That is, any
        attributes that are containers are shared between the new graph
        and the original. This is exactly what `dict.copy()` provides.
        You can obtain this style copy using:

            >>> G = nx.path_graph(5)
            >>> H = G.copy()
            >>> H = G.copy(as_view=False)
            >>> H = nx.Graph(G)
            >>> H = G.fresh_copy().__class__(G)

        Fresh Data -- For fresh data, the graph structure is copied while
        new empty data attribute dicts are created. The resulting graph
        is independent of the original and it has no edge, node or graph
        attributes. Fresh copies are not enabled. Instead use:

            >>> H = G.fresh_copy()
            >>> H.add_nodes_from(G)
            >>> H.add_edges_from(G.edges)

        View -- Inspired by dict-views, graph-views act like read-only
        versions of the original graph, providing a copy of the original
        structure without requiring any memory for copying the information.

        See the Python copy module for more information on shallow
        and deep copies, https://docs.python.org/2/library/copy.html.

        Parameters
        ----------
        as_view : bool, optional (default=False)
            If True, the returned graph-view provides a read-only view
            of the original graph without actually copying any data.

        Returns
        -------
        G : Graph
            A copy of the graph.

        See Also
        --------
        to_directed: return a directed copy of the graph.

        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> H = G.copy()

        """
        raise NotImplementedError
        # if as_view is True:
        #     return nx.graphviews.GraphView(self)
        # G = self.fresh_copy()
        # G.graph.update(self.graph)
        # G.add_nodes_from((n, d.copy()) for n, d in self._node.items())
        # G.add_edges_from((u, v, datadict.copy())
        #                  for u, nbrs in self._adj.items()
        #                  for v, datadict in nbrs.items())
        # return G

    def to_directed(self, as_view=False):
        """Return a directed representation of the graph.

        Returns
        -------
        G : DiGraph
            A directed graph with the same name, same nodes, and with
            each edge (u, v, data) replaced by two directed edges
            (u, v, data) and (v, u, data).

        See Also
        -----
        Graph.to_directed
        """
        raise nx.NetworkXError("Map doesn't support to_directed")
        # if as_view is True:
        #     return nx.graphviews.DiGraphView(self)
        # # deepcopy when not a view
        # from networkx import DiGraph
        # G = DiGraph()
        # G.graph.update(deepcopy(self.graph))
        # G.add_nodes_from((n, deepcopy(d)) for n, d in self._node.items())
        # G.add_edges_from((u, v, deepcopy(data))
        #                  for u, nbrs in self._adj.items()
        #                  for v, data in nbrs.items())
        # return G

    def to_undirected(self, as_view=False):
        """Return an undirected copy of the graph.

        Parameters
        ----------
        as_view : bool (optional, default=False)
          If True return a view of the original undirected graph.

        Returns
        -------
        G : Graph/MultiGraph
            A deepcopy of the graph.

        See Also
        --------
        Graph, copy, add_edge, add_edges_from, Graph.to_undirected
        """
        raise nx.NetworkXError("Map doesn't support to_undirected")
        # if as_view is True:
        #     return nx.graphviews.GraphView(self)
        # # deepcopy when not a view
        # G = Graph()
        # G.graph.update(deepcopy(self.graph))
        # G.add_nodes_from((n, deepcopy(d)) for n, d in self._node.items())
        # G.add_edges_from((u, v, deepcopy(d))
        #                  for u, nbrs in self._adj.items()
        #                  for v, d in nbrs.items())
        # return G