import networkx as nx
from nose.tools import assert_equals
from nose.tools import ok_
from networkx.algorithms.approximation import treewidth_min_degree


def is_tree_decomp(graph, decomp):
    isTreeDecomp = True
    for x in graph.nodes():
        appearOnce = False
        for bag in decomp.nodes():
            if (x in bag):
                appearOnce = True
        ok_(appearOnce)

    # Check if each connected pair of nodes appears at least once together in a bag
    for (x, y) in graph.edges():
        appearTogether = False
        for bag in decomp.nodes():
            if (x in bag and y in bag):
                appearTogether = True
        ok_(appearTogether)

        # Check if the nodes associated with vertex v form a connected subset of T
    for v in graph.nodes():
        subset = []
        for bag in decomp.nodes():
            if (v in bag):
                subset.append(bag)
        subGraph = decomp.subgraph(subset)
        ok_(nx.is_connected(subGraph))


class TestTreewidth(object):
    def test_tree_decomposition(self):
        # Test if the returned decomposition is a valid decomposition for a selection of various graphs
        G = nx.petersen_graph()
        treewidth, decomp = treewidth_min_degree(G)
        is_tree_decomp(G, decomp)
        # Check if each vertex appears at least once in a bag

    def test_min_degreee(self):
        G = nx.Graph()
        # The order of removal should be (with [] denoting any order of the containing nodes) [1,2,4]3[5,6,7], resulting in treewdith 2 for the heuristic
        G.add_edge(1, 3)
        G.add_edge(4, 3)
        G.add_edge(2, 3)
        G.add_edge(3, 5)
        G.add_edge(5, 6)
        G.add_edge(5, 7)
        G.add_edge(6, 7)
        treewidth, decomp = treewidth_min_degree(G)
        assert_equals(treewidth, 2)
