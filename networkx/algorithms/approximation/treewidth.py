# -*- coding: utf-8 -*-
import sys

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


# version without recursion and copying, rest is still the same as in tree_decomp
def tree_decomp_unopt(G, heuristic):
    # Copy graph because the algorithm modifies it
    G = G.copy()
    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []
    timings = []
    elim_node = heuristic(G)
    while elim_node is not None:

        # Connect all neighbours with each other
        neighbors = set(G.neighbors(elim_node))
        for n in neighbors:
            for m in neighbors:
                if (n != m):
                    G.add_edge(n, m)

        # remove node from graph and push on stack (including its neighbors)
        G.remove_node(elim_node)
        node_stack.append((elim_node, neighbors))
        timings.append(time.time())
        timings.append(time.time())

        # get next node to be removed according to heuristic
        elim_node = heuristic(G)

    timings.append(time.time())
    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    decomp.add_node(frozenset(G.nodes))

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

        # Add edge to decomposition (implicitly also adds the new node)
        decomp.add_edge(old_bag, neighbors)
    timings.append(time.time())
    return decomp, timings


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
