# -*- coding: utf-8 -*-
import sys,timeit

import numpy as np

import networkx as nx
import time

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]



# Returns a tuple: (treewidth: int, decomposition: Graph)
def treewidth_min_degree(G):
    return treewidth_decomp(G, min_degree_heuristic)


# Returns a tuple: (treewidth: int, decomposition: Graph)
def treewidth_min_fill_in(G):
    return treewidth_decomp(G, min_fill_in_heuristic)


# Returns the node that needs to be removed next or None (if the abort condition is met)
def min_degree_heuristic(G):
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


# Returns the node that needs to be removed next or None (if the abort condition is met)
def min_degree_heuristic_numpy(adj_array):
    degree_array = np.count_nonzero(adj_array,0)
    nr_non_deleted = np.count_nonzero(degree_array)
    if nr_non_deleted == 0:
        # All nodes removed
        return None

    min_degree = np.amin(degree_array[degree_array>0])

    if min_degree == nr_non_deleted: #Note min_degree always contains self loops
        return None
    else:
        return np.where(degree_array==min_degree)[0][0]


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
    # Copy graph because the algorithm modifies it
    G = G.copy()
    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    # Maps a node n to a lists of nodes l
    # Later when bag(n) is created every bag(node in l) might need to connect to bag(n)
    old_bag_notify = {}

    elim_node = heuristic(G)
    while elim_node is not None:

        # Connect all neighbours with each other
        neighbors = list(G[elim_node])
        for i, node_1 in enumerate(neighbors[:-1]):  # Iterate neighbors excluding last element
            for node_2 in neighbors[i+1:]:  # Iterate neighbors after node_1
                if not G.has_edge(node_1, node_2):
                    G.add_edge(node_1, node_2)

        # Some neighbor will create the bag that elim_node later needs to connect to
        for node in neighbors:
            if node in old_bag_notify:
                old_bag_notify[node].append(elim_node)
            else:
                old_bag_notify[node] = [elim_node]

        # remove node from graph and push on stack (including its neighbors)
        G.remove_node(elim_node)
        node_stack.append((elim_node, neighbors))

        # get next node to be removed according to heuristic
        elim_node = heuristic(G)

    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    first_bag = frozenset(G.nodes)
    decomp.add_node(first_bag)

    treewidth = len(first_bag) - 1

    # Maps node n to the bag to which bag(n) needs to connect
    old_bag_connection = {}
    for node in G.nodes:
        old_bag_connection[node] = first_bag

    while node_stack:
        # get node and its neighbors from the stack
        (curr_node, neighbors) = node_stack.pop()

        # find a bag the neighbors are in (or default to the first created bag)
        if curr_node in old_bag_connection:
            old_bag = old_bag_connection[curr_node]
        else:
            old_bag = first_bag

        # Create new node for decomposition
        neighbors.append(curr_node)
        new_bag = frozenset(neighbors)
        # Update treewidth
        treewidth = max(treewidth, len(new_bag)-1)

        # If this node was the first in a created clique to get deleted the created bag is the old_bag from this node
        if curr_node in old_bag_notify:
            for old_neighbor_node in old_bag_notify[curr_node]:
                # set (possibly override) the bag of old_neighbor_node should connect to
                old_bag_connection[old_neighbor_node] = new_bag

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, new_bag)
    return treewidth, decomp


# Calculates tree width decomposition using the passed heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp_numpy(G: nx.Graph, heuristic):
    adj_matrix = nx.to_numpy_array(G,dtype=np.bool_)
    N = len(adj_matrix)
    #Add self loops (so all nodes have at least degree 1, and removed nodes have degree 0)
    adj_matrix[range(N), range(N)] = 1

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    # Maps a node n to a lists of nodes l
    # Later when bag(n) is created every bag(node in l) might need to connect to bag(n)
    old_bag_notify = {}

    elim_node = heuristic(adj_matrix)
    while elim_node is not None:

        # Connect all neighbours with each other
        neighbor_mask = adj_matrix[elim_node]
        neighbors = np.where(neighbor_mask)[0]
        neighbors_x = np.tile(neighbors,len(neighbors))
        neighbors_y = np.repeat(neighbors, len(neighbors))
        # Add new edges
        adj_matrix[neighbors_x,neighbors_y] = 1


        # Some neighbor will create the bag that elim_node later needs to connect to
        for node in neighbors:
            if node != elim_node:
                if node in old_bag_notify:
                    old_bag_notify[node].append(elim_node)
                else:
                    old_bag_notify[node] = [elim_node]

        # remove node from graph and push on stack (including its neighbors)
        adj_matrix[elim_node]   = 0
        adj_matrix[:,elim_node] = 0
        node_stack.append((elim_node, neighbors))

        # get next node to be removed according to heuristic
        elim_node = heuristic(adj_matrix)

    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    first_bag = frozenset(adj_matrix.nonzero()[0])
    decomp.add_node(first_bag)

    treewidth = len(first_bag) - 1

    # Maps node n to the bag to which bag(n) needs to connect
    old_bag_connection = {}
    for node in G.nodes:
        old_bag_connection[node] = first_bag

    while node_stack:
        # get node and its neighbors from the stack
        (curr_node, neighbors) = node_stack.pop()

        # find a bag the neighbors are in (or default to the first created bag)
        if curr_node in old_bag_connection:
            old_bag = old_bag_connection[curr_node]
        else:
            old_bag = first_bag

        # Create new node for decomposition
        neighbors = list(neighbors)
        neighbors.append(curr_node)
        new_bag = frozenset(neighbors)
        # Update treewidth
        treewidth = max(treewidth, len(new_bag)-1)

        # If this node was the first in a created clique to get deleted the created bag is the old_bag from this node
        if curr_node in old_bag_notify:
            for old_neighbor_node in old_bag_notify[curr_node]:
                # set (possibly override) the bag of old_neighbor_node should connect to
                old_bag_connection[old_neighbor_node] = new_bag

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, new_bag)
    return treewidth, decomp


if __name__ == '__main__':
    G = nx.fast_gnp_random_graph(1500,0.001)
    NUM_EXEC = 10
    print(timeit.timeit('tw(G,heur)', globals={'G': G, 'tw': treewidth_decomp_numpy, 'heur': min_degree_heuristic_numpy},number=NUM_EXEC))
    print(timeit.timeit('tw(G,heur)',globals={'G':G,'tw':treewidth_decomp,'heur':min_degree_heuristic},number=NUM_EXEC))
