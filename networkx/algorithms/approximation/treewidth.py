# -*- coding: utf-8 -*-
"""Functions for computing treewidth decomposition.
   
   Treewidth of an undirected graph is a number associated with the graph. It can be defined as the size of the largest vertex set (bag) in a tree decomposition of the graph minus one. 

`Wikipedia: Treewidth <https://en.wikipedia.org/wiki/Treewidth>`_

   The notions of treewidth and tree decomposition have gained their attractiveness partly because many graph and network problems that are intractable (e.g., NP-hard) on arbitrary graphs become efficiently solvable (e.g., with a linear time algorithm) when the treewidth of the input graphs is bounded by a constant.

    Hans L. Bodlaender and Arie M. C. A. Koster. 2010. "Treewidth computations I.Upper bounds". Inf. Comput. 208, 3 (March 2010),259-275. DOI=http://dx.doi.org/10.1016/j.ic.2009.03.008 

    Hand L. Bodlaender. "Discovering Treewidth". institute of information and computing sciences, utrecht university. technical report UU-CS-2005-018. www.cs.uu.nl

"""
import sys

import networkx as nx
from networkx.utils import not_implemented_for
from heapq import heappush, heappop, heapify

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]


# Returns a tuple: (treewidth: int, decomposition: Graph)
@not_implemented_for('directed')
@not_implemented_for('multigraph')
def treewidth_min_degree(G):
    """ Returns a treewidth decomposition using the Minimum Degree heuristic.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
          2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """
    return treewidth_decomp(G, MinDegreeHeuristic)


# Returns a tuple: (treewidth: int, decomposition: Graph)
@not_implemented_for('directed')
@not_implemented_for('multigraph')
def treewidth_min_fill_in(G):
    """ Returns a treewidth decomposition using the Minimum Fill-in heuristic.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
        2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """
    return treewidth_decomp(G, MinFillInHeuristic)


class MinDegreeHeuristic:
    # TODO: Update documentation
    def __init__(self, graph):
        self._graph = graph

        # A collection of nodes that have to be updated in the heap before each iteration
        self._update_nodes = []

        # self._degreeq is heapq with 2-tuples (degree,node)
        self._degreeq = []
        # build heap with initial degrees
        for n in graph:
            self._degreeq.append((len(graph[n]), n))
        heapify(self._degreeq)

    def __iter__(self):
        return self

    # Implement next method for backwards compatibility with python 2
    def next(self):
        return self.__next__()

    def __next__(self):
        # Update nodes in self._update_nodes
        for n in self._update_nodes:
            # insert changed degrees into degreeq
            heappush(self._degreeq, (len(self._graph[n]), n))

        while self._degreeq:
            # get the next (minimum degree) node
            (min_degree, elim_node) = heappop(self._degreeq)
            if not elim_node in self._graph or len(self._graph[elim_node]) != min_degree:
                # Outdated entry in degreeq
                continue
            elif min_degree == len(self._graph) - 1:
                # Fully connected: Abort condition
                raise StopIteration

            # Remember to update nodes in the heap before getting the next node
            self._update_nodes = self._graph[elim_node]
            return elim_node

        # The heap is empty: Abort
        raise StopIteration


class MinFillInHeuristic:
    """Returns the node from the graph, where the number of edges added  when
    turning the neighbourhood of the chosen node into clique is small as possible.
    Parameters
    ----------
    G : Graph
    Returns
    -------
    min_fill_node : string, integers or hashable Python object (except None)
        The node from the graph, for which, when it is deleted from the graph and
        its neighbourhood is turned into clique, the number of edges added is
        small as possible.
    Notes
    -----
    This algorithm computes the node with 'min fill in' in the graph 'G'.
    The running time of the algorithm is O(V*V*V) and it uses constant
    additional memory.
    References
    ----------
    .. [1] K. Wang, Z. Lu, and J. Hicks
           *Treewidth*.
           http://web.eecs.utk.edu/~cphillip/cs594_spring2015_projects/treewidth.pdf
    """
    #TODO: Update documentation
    def __init__(self, graph):
        self._graph = graph

    def __iter__(self):
        return self

    # Implement next method for backwards compatibility with python 2
    def next(self):
        return self.__next__()

    def __next__(self):
        if len(self._graph) == 0:
            raise StopIteration

        min_fill_in_node = None

        min_fill_in = sys.maxsize

        # create sorted list of (degree, node)
        degree_list = [(len(self._graph[node]), node) for node in self._graph]
        degree_list.sort()

        # abort condition
        min_degree = degree_list[0][0]
        if min_degree == len(self._graph) - 1:
            raise StopIteration

        for (degree, node) in degree_list:
            num_fill_in = 0
            # Convert to list in order to access by index
            neighbors = list(self._graph[node])
            for i in range(len(neighbors) - 1):
                for j in range(i + 1, len(neighbors)):
                    if not neighbors[j] in self._graph[neighbors[i]]:
                        num_fill_in += 1
                        # prune if this can't be min-fill-in node anymore
                        if num_fill_in >= min_fill_in:
                            break
                else:
                    continue  # executed if no break
                break

            if num_fill_in < min_fill_in:
                if num_fill_in == 0:
                    return node
                min_fill_in = num_fill_in
                min_fill_in_node = node

        return min_fill_in_node


# Calculates tree width decomposition using the passed heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp(G, heuristic_class):
    """Returns a treewidth decomposition using the passed heuristic.

    Parameters
    ----------
    G : NetworkX graph
    heuristic_class : iterator class

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
        2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """
    # TODO: Update documentation regarding heuristic_class
    # make dict-of-sets structure
    graph = {}
    for u in G:
        graph[u] = set()
        for v in G[u]:
            if u != v:  # ignore self-loop
                graph[u].add(v)

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    # instanciate a heuristic_iterator
    heuristic_iterator = heuristic_class(graph)

    for elim_node in heuristic_iterator:
        # Connect all neighbours with each other
        neighbors = graph[elim_node]
        for n in neighbors:
            for m in neighbors:
                if n != m and not m in graph[n]:
                    graph[n].add(m)
                    graph[m].add(n)

        # push node and its current neighbors on stack
        node_stack.append((elim_node, neighbors))

        # remove node from graph
        for u in graph:
            if elim_node in graph[u]:
                graph[u].remove(elim_node)
        del graph[elim_node]

    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    first_bag = frozenset(graph.keys())
    decomp.add_node(first_bag)

    treewidth = len(first_bag) - 1

    while node_stack:
        # get node and its neighbors from the stack
        (curr_node, neighbors) = node_stack.pop()

        # find a bag the neighbors are in
        old_bag = None
        for bag in decomp.nodes:
            if neighbors <= bag:
                old_bag = bag
                break
        if old_bag == None:
            # no old_bag was found: just connect to the first_bag
            old_bag = first_bag

        # Create new node for decomposition
        neighbors.add(curr_node)
        new_bag = frozenset(neighbors)

        # Update treewidth
        treewidth = max(treewidth, len(new_bag) - 1)

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, new_bag)

    print(treewidth, decomp.edges)
    return treewidth, decomp
