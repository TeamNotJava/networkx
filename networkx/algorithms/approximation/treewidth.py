# -*- coding: utf-8 -*-
"""Functions for computing treewidth decomposition.
   
Treewidth of an undirected graph is a number associated with the graph.
It can be defined as the size of the largest vertex set (bag) in a tree
decomposition of the graph minus one.

`Wikipedia: Treewidth <https://en.wikipedia.org/wiki/Treewidth>`_

The notions of treewidth and tree decomposition have gained their
attractiveness partly because many graph and network problems that are
intractable (e.g., NP-hard) on arbitrary graphs become efficiently
solvable (e.g., with a linear time algorithm) when the treewidth of the
input graphs is bounded by a constant [1]_ [2]_.

There are two classes which contain implementations of different heuristics for
computing tree decomposition: :class:`MinDegreeHeuristic` and 
:class:` min_fill_in_heuristic`.
   
:class:`MinDegreeHeuristic`
    Returns a treewidth decomposition using the Minimum Degree heuristic.
    The heuristic chooses the nodes according to their degree
    (number of neighours), i.e., first the node with the lowest degree is
    chosen, then the graph is updated and the correspondig node is
    removed. Next, a new node with the lowest degree is chosen,
    and so on.
   
        
:class:` min_fill_in_heuristic`
    Returns the node from the graph, where the number of edges added  when
    turning the neighbourhood of the chosen node into clique is as small as
    possible. This algorithm chooses the nodes using the Minimum Fill-In
    heuristic. The running time of the algorithm is :math:`O(V^3)` and it uses
    additional constant memory [3]_.
       
   
.. [1] Hans L. Bodlaender and Arie M. C. A. Koster. 2010. "Treewidth computations
      I.Upper bounds". Inf. Comput. 208, 3 (March 2010),259-275.
      http://dx.doi.org/10.1016/j.ic.2009.03.008

.. [2] Hand L. Bodlaender. "Discovering Treewidth". Institute of Information and
      Computing Sciences, Utrecht University. Technical Report UU-CS-2005-018.
      http://www.cs.uu.nl
   
.. [3] K. Wang, Z. Lu, and J. Hicks *Treewidth*.
      http://web.eecs.utk.edu/~cphillip/cs594_spring2015_projects/treewidth.pdf

"""

import sys

import networkx as nx
from networkx.utils import not_implemented_for
from heapq import heappush, heappop, heapify
import itertools

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def treewidth_min_degree(G):
    """ Returns a treewidth decomposition using the Minimum Degree heuristic. The
        heuristic chooses the nodes according to their degree, i.e., first the node
        with the lowest degree is chosen, then the graph is updated and the correspondig
        node is removed. Next, a new node with the lowest degree is chosen, and so on.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
          2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """

    deg_heuristic = MinDegreeHeuristic(G)
    return treewidth_decomp(G, lambda graph: deg_heuristic.best_node(graph))


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def treewidth_min_fill_in(G):
    """ Returns a treewidth decomposition using the Minimum Fill-in heuristic. The
    heuristic chooses a node from the graph, where the number of edges added when
    turning the neighbourhood of the chosen node into clique is small as possible.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
        2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """
    return treewidth_decomp(G,  min_fill_in_heuristic)


class MinDegreeHeuristic:
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

    def best_node(self,graph):
         # Update nodes in self._update_nodes
        for n in self._update_nodes:
            # insert changed degrees into degreeq
            heappush(self._degreeq, (len(graph[n]), n))
        while self._degreeq:
            # get the next (minimum degree) node
            (min_degree, elim_node) = heappop(self._degreeq)
            if elim_node not in graph or len(graph[elim_node]) != min_degree:
                # Outdated entry in degreeq
                continue
            elif min_degree == len(graph) - 1:
                # Fully connected: Abort condition
                return None

            # Remember to update nodes in the heap before getting the next node
            self._update_nodes = graph[elim_node]
            return elim_node

        # The heap is empty: Abort
        return None




def min_fill_in_heuristic(graph):
    if len(graph) == 0:
        return None

    min_fill_in_node = None

    min_fill_in = sys.maxsize

    # create sorted list of (degree, node)
    degree_list = [(len(graph[node]), node) for node in graph]
    degree_list.sort()

    # abort condition
    min_degree = degree_list[0][0]
    if min_degree == len(graph) - 1:
        return None

    for (_, node) in degree_list:
        num_fill_in = 0

        nbrs = graph[node]
        for nbr in nbrs:
            # count how many nodes in nbrs current nbr is not connected to
            # - 1 for the node itself
            num_fill_in += len(nbrs - graph[nbr]) - 1
            if num_fill_in >= 2 * min_fill_in:
                break

        # divide by 2 because of double counting
        num_fill_in /= 2

        if num_fill_in < min_fill_in: # Update min-fill-in node
            if num_fill_in == 0:
                return node
            min_fill_in = num_fill_in
            min_fill_in_node = node

    return min_fill_in_node


def treewidth_decomp(G, heuristic=min_fill_in_heuristic):
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
    
    # make dict-of-sets structure

    graph = {n:set(G[n]) - set([n]) for n in G}

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    # instantiate a heuristic_iterator
    #heuristic_iterator = heuristic_class(graph)
    elim_node=heuristic(graph)
    while not elim_node is None:
        # Connect all neighbours with each other
        nbrs = graph[elim_node]
        for u, v in itertools.permutations(nbrs, 2):
            if v not in graph[u]:
                graph[u].add(v)

        # push node and its current neighbors on stack
        node_stack.append((elim_node, nbrs))

        # remove node from graph
        for u in graph:
            if elim_node in graph[u]:
                graph[u].remove(elim_node)
        del graph[elim_node]
        elim_node=heuristic(graph)
    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    first_bag = frozenset(graph.keys())
    decomp.add_node(first_bag)

    treewidth = len(first_bag) - 1

    while node_stack:
        # get node and its neighbors from the stack
        (curr_node, nbrs) = node_stack.pop()

        # find a bag the neighbors are in
        old_bag = None
        for bag in decomp.nodes:
            if nbrs <= bag:
                old_bag = bag
                break
        if old_bag == None:
            # no old_bag was found: just connect to the first_bag
            old_bag = first_bag

        # Create new node for decomposition
        nbrs.add(curr_node)
        new_bag = frozenset(nbrs)

        # Update treewidth
        treewidth = max(treewidth, len(new_bag) - 1)

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, new_bag)

    return treewidth, decomp
