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

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]




# Returns a tuple: (treewidth: int, decomposition: Graph)
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
    return treewidth_decomp(G, min_degree_heuristic)


# Returns a tuple: (treewidth: int, decomposition: Graph)
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
    return treewidth_decomp(G, min_fill_in_heuristic)


# Returns the node that needs to be removed next or None (if the abort condition is met)
def min_degree_heuristic(G):
    """ Returns the node that needs to be removed next or None according to the Minimum Degree heuristic that selects the vertex with the minimum number of unselected neighbours.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Node : NetworkX graph
         
    References
    ----------
    .. [2] Hand L. Bodlaender. "Discovering Treewidth". institute of information and computing sciences, utrecht university. technical report UU-CS-2005-018. www.cs.uu.nl
    """
    min_degree = sys.maxsize
    min_node = None
    for (node, degree) in G.degree:
        if degree < min_degree:
            if degree <= 1:
                # Return early
                return node
            min_node = node
            min_degree = degree

    if min_degree == G.number_of_nodes() - 1:
        # Fully connected: Abort condition
        return None
    else:
        return min_node


# Returns the node with minimum degree or None (if the abort condition is met)
def min_fill_in_heuristic(G):
    """Returns the node with minimum degree or None according to the Minimum Fill-in heuristic that selects a vertex which geves the minimum number of added fill-in edges for the current step.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    Node : NetworkX graph

    References
    ----------
    .. [2] Hand L. Bodlaender. "Discovering Treewidth". institute of information and computing sciences, utrecht university. technical report UU-CS-2005-018. www.cs.uu.nl
    """
    candidate_node = None
    # Still keep track of min_degree to abort earlier
    min_degree = sys.maxsize
    min_fill_in = sys.maxsize

    for (node, degree) in G.degree:
        min_degree = min(min_degree, degree)
        num_fill_in = 0
        # Convert to list in order to access by index
        neighbors = list(G.neighbors(node))
        for i in range(len(neighbors) - 1):
            for j in range(i + 1, len(neighbors)):
                if not G.has_edge(neighbors[i], neighbors[j]):
                    num_fill_in += 1

        if num_fill_in < min_fill_in:
            if num_fill_in == 0:
                return node
            min_fill_in = num_fill_in
            candidate_node = node

    if min_degree == G.number_of_nodes() - 1:
        # Fully connected: Abort condition
        return None
    else:
        return candidate_node


# Calculates tree width decomposition using the passed heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp(G, heuristic):
    """Returns a treewidth decomposition using the passed heuristic.

    Parameters
    ----------
    G : NetworkX graph
    heuristic : function

    Returns
    -------
    Treewidth decomposition : (int, Graph) tuple
        2-tuple with treewidth and the corresponding decomposed tree (NetworkX graph).
    """
    # copy so the original graph is not modified
    G = G.copy()

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    elim_node = heuristic(G)
    while elim_node is not None:

        # Connect all neighbours with each other
        neighbors = set(G.neighbors(elim_node))
        for n in neighbors:
            for m in neighbors:
                if n != m and not G.has_edge(n, m):
                    G.add_edge(n, m)

        # remove node from graph and push on stack (including its neighbors)
        G.remove_node(elim_node)
        node_stack.append((elim_node, neighbors))

        # get next node to be removed according to heuristic
        elim_node = heuristic(G)

    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    new_bag = frozenset(G.nodes)
    decomp.add_node(new_bag)

    current_treewidth = len(new_bag) - 1

    while node_stack:
        # get node and its neighbors from the stack
        (curr_node, neighbors) = node_stack.pop()
        # find a bag the neighbors are in
        old_bag = None
        for bag in decomp.nodes:
            if neighbors <= bag:
                old_bag = bag
                break

        # Create new node for decomposition
        neighbors.add(curr_node)
        neighbors = frozenset(neighbors)

        # Update treewidth
        current_treewidth = max(current_treewidth, len(neighbors)-1)

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, neighbors)

    return current_treewidth, decomp
