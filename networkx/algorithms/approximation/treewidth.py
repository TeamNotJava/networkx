# -*- coding: utf-8 -*-
import sys

from heapq import heappush, heappop, heapify
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


# Calculates tree width decomposition using the passed heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp(G, heuristic):
    G = G.copy()

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    elim_node = heuristic(G)
    while elim_node is not None:
        # Connect all neighbours with each other
        neighbors = set(G.neighbors(elim_node))
        for n in neighbors:
            for m in neighbors:
                if (n != m) and not G.has_edge(n, m):
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

# Calculates tree width decomposition according to the minimum degree heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp_min_degree(G):
    G = G.copy()

    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    push = heappush
    pop = heappop
    # degreeq is heapq with 2-tuples (degree,node)
    degreeq = []

    # build heap with initial degrees
    for (n, degree) in G.degree:
        degreeq.append((degree, n))
    heapify(degreeq)

    while degreeq:
        # get the next (minimum degree) node
        (min_degree, elim_node) = pop(degreeq)
        if not G.has_node(elim_node) or G.degree[elim_node] != min_degree:
            # Outdated entry in degreeq
            continue
        elif min_degree == G.number_of_nodes() - 1:
            # Fully connected: Abort condition
            break

        # Connect all neighbours with each other
        neighbors = set(G.neighbors(elim_node))
        changed_degree = set(G.neighbors(elim_node))
        for n in neighbors:
            for m in neighbors:
                if (n != m) and not G.has_edge(n, m):
                    G.add_edge(n, m)
                    changed_degree.add(n)
                    changed_degree.add(m)

        # remove node from graph and push on stack (including its neighbors)
        G.remove_node(elim_node)
        node_stack.append((elim_node, neighbors))

        # insert changed degrees into degreeq
        for n in changed_degree:
            push(degreeq, (G.degree[n], n))

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


# Calculates tree width decomposition according to the minimum degree heuristic
# Returns tuple: (treewidth: int, decomposition: Graph)
def treewidth_decomp_min_degree_faster(G):
    #G = G.copy()
    graph = {}
    for u in G:
        graph[u] = set()
        for v in G[u]:
            graph[u].add(v)


    # stack where nodes and their neighbors are pushed in the order they are selected by the heuristic
    node_stack = []

    push = heappush
    pop = heappop
    # degreeq is heapq with 2-tuples (degree,node)
    degreeq = []

    # build heap with initial degrees
    #for (n, degree) in G.degree:
    #    degreeq.append((degree, n))
    for n in graph:
        degreeq.append((len(graph[n]), n))
    heapify(degreeq)

    while degreeq:
        # get the next (minimum degree) node
        (min_degree, elim_node) = pop(degreeq)
        #if not G.has_node(elim_node) or G.degree[elim_node] != min_degree:
        if not elim_node in graph or len(graph[elim_node]) != min_degree:
            # Outdated entry in degreeq
            continue
        #elif min_degree == G.number_of_nodes() - 1:
        elif min_degree == len(graph) - 1:
            # Fully connected: Abort condition
            break

        # Connect all neighbours with each other
        #neighbors = set(G.neighbors(elim_node))
        neighbors = graph[elim_node]
        changed_degree = neighbors
        for n in neighbors:
            for m in neighbors:
                #if (n != m) and not G.has_edge(n, m):
                if n != m and not m in graph[n]:
                    #G.add_edge(n, m)
                    graph[n].add(m)
                    graph[m].add(n)
                    changed_degree.add(n)
                    changed_degree.add(m)

        # remove node from graph and push on stack (including its neighbors)
        #G.remove_node(elim_node)
        for u in graph:
            if elim_node in graph[u]:
                graph[u].remove(elim_node)
        graph.pop(elim_node, None)
        node_stack.append((elim_node, neighbors))

        # insert changed degrees into degreeq
        for n in changed_degree:
            #push(degreeq, (G.degree[n], n))
            push(degreeq, (len(graph[n]), n))

    # The abort condition is met. Put all nodes into one bag.
    decomp = nx.Graph()
    #new_bag = frozenset(G.nodes)
    new_bag = frozenset(graph.keys())
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
