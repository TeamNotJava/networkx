# -*- coding: utf-8 -*-
import sys

from networkx import Graph

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]



# Returns a tuple: (treewidth: int, decomposition: Graph)
def treewidth_min_degree(G):
    return treewidth_decomp(G, min_degree_node_heuristic)


# Returns a tuple: (treewidth: int, decomposition: Graph)
def treewidth_min_fill_in(G):
    return treewidth_decomp(G, min_fill_in_heuristic)


def min_degree_node_heuristic(G: Graph):
    """Returns the node from the graph with minimum degree.

        Parameters
        ----------
        G : Graph

        Returns
        -------
        min_degree_node : string, integers or hashable Python object (except None)
            The node from the graph with lowest degree. If the graph is fully connected
            then the function returns None.

        Notes
        -----
        This algorithm computes the node with lowest number of edges in the graph 'G'.
        The running time of the algorithm is O(V) and it uses constant
        additional memory.

        References
        ----------
        .. [1] K. Wang, Z. Lu, and J. Hicks
               *Treewidth*.
               http://web.eecs.utk.edu/~cphillip/cs594_spring2015_projects/treewidth.pdf

        """

    min_degree = sys.maxsize
    min_degree_node = None
    for (node, degree) in G.degree:
        if degree < min_degree:
            if degree <= 1:
                # Return early
                return node
            min_degree_node = node
            min_degree = degree

    if min_degree == G.number_of_nodes() - 1:
        # Fully connected: Abort condition
        return None
    else:
        return min_degree_node



def min_fill_in_heuristic_faster(G):
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

    if len(G) == 0:
        return None

    min_fill_in_node = None

    min_fill_in = sys.maxsize

    # create sorted list of (degree, node)
    degree_list = [(degree, node) for (node, degree) in G.degree]
    degree_list.sort()

    # abort condition
    min_degree = degree_list[0][0]
    if min_degree == len(G) - 1:
        return None

    for (degree, node) in degree_list:
        min_degree = min(min_degree, degree)
        num_fill_in = 0
        # Convert to list in order to access by index
        neighbors = list(G.neighbors(node))
        for i in range(len(neighbors) - 1):
            for j in range(i + 1, len(neighbors)):
                if not G.has_edge(neighbors[i], neighbors[j]):
                    num_fill_in += 1
                    # prune if this can't be min-fill-in node anymore
                    if num_fill_in >= min_fill_in:
                        break
            else:
                continue # executed if no break
            break

        if num_fill_in < min_fill_in:
            if num_fill_in == 0:
                return node
            min_fill_in = num_fill_in
            min_fill_in_node = node

    return min_fill_in_node

# Returns the node with minimum degree or None (if the abort condition is met)
def min_fill_in_heuristic(G):
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
    decomp = Graph()
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
