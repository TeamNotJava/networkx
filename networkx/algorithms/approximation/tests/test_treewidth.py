#-*- coding: utf-8 -*-
#    Copyright (C) 2018 by
#    Rudolf-Andreas Floren <rudi.floren@gmail.com>
#    Dominik Meier <dominik.meier@rwth-aachen.de>
#    All rights reserved.
#    BSD license.
import networkx as nx
from nose.tools import assert_equals
from nose.tools import ok_
from networkx.algorithms.approximation import treewidth_min_degree
from networkx.algorithms.approximation import treewidth_min_fill_in
"""Unit tests for the :mod:`networkx.algorithms.approximation.treewidth`
module.

"""


def is_tree_decomp(graph, decomp):
    """Check if the given decomposition tree is a decomposition tree of the given graph"""
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


class TestTreewidthMinDegree(object):
    """Unit tests for the :mod:`networkx.algorithms.approximation.treewidth_min_degree`
    function.

    """

    def setUp(self):
        """Setup for different kinds of trees"""
        self.smallTree = nx.Graph()
        self.smallTree.add_edge(1, 3)
        self.smallTree.add_edge(4, 3)
        self.smallTree.add_edge(2, 3)
        self.smallTree.add_edge(3, 5)
        self.smallTree.add_edge(5, 6)
        self.smallTree.add_edge(5, 7)
        self.smallTree.add_edge(6, 7)

    def test_tree_decomposition(self):
        """Test if the returned decomposition is a valid decomposition for a selection 
        of various graphs
        """
        G = nx.petersen_graph()
        _, decomp = treewidth_min_degree(G)
        is_tree_decomp(G, decomp)
        # Check if each vertex appears at least once in a bag

    def test_small_tree_treewidth(self):
        """Test if the computed treewidth of the known self.smallTree is 2.
        As we know which value we can expect from our heuristic, values other than two are regressions

        """
        G = self.smallTree
        # The order of removal should be (with [] denoting any order of the containing nodes) [1,2,4]3[5,6,7], resulting in treewdith 2 for the heuristic
        
        treewidth, _ = treewidth_min_fill_in(G)
        assert_equals(treewidth, 2)


    def test_mid_tree_treewidth(self):
        pass

    def test_large_tree_treewidth(self):
        pass


class TestTreewidthMinFillIn(object):
    """Unit tests for the :mod:`networkx.algorithms.approximation.treewidth_min_fill_in`
    function.

    """

    def setUp(self):
        """Setup for different kinds of trees"""
        self.smallTree = nx.Graph()
        self.smallTree.add_edge(1, 3)
        self.smallTree.add_edge(4, 3)
        self.smallTree.add_edge(2, 3)
        self.smallTree.add_edge(3, 5)
        self.smallTree.add_edge(5, 6)
        self.smallTree.add_edge(5, 7)
        self.smallTree.add_edge(6, 7)

    def test_tree_decomposition(self):
        """Test if the returned decomposition is a valid decomposition for a selection 
        of various graphs
        """
        G = nx.petersen_graph()
        _, decomp = treewidth_min_fill_in(G)
        is_tree_decomp(G, decomp)
        # Check if each vertex appears at least once in a bag

    def test_small_tree_treewidth(self):
        """Test if the computed treewidth of the known self.smallTree is 2"""
        G = self.smallTree
        # The order of removal should be (with [] denoting any order of the containing nodes) [1,2,4]3[5,6,7], resulting in treewdith 2 for the heuristic
        
        treewidth, _ = treewidth_min_degree(G)
        assert_equals(treewidth, 2)

    def test_mid_tree_treewidth(self):
        pass

    def test_large_tree_treewidth(self):
        pass