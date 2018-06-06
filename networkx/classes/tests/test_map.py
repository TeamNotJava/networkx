from nose.tools import assert_equal
from nose.tools import assert_is
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises

import networkx as nx
from networkx.testing.utils import *
from test_graph import BaseGraphTester

class TestMap(object):
    def setup(self):
        self.w_3 = nx.Map()
        self.w_3.add_edge(0, 1)
        self.w_3.add_edge(1, 2)
        self.w_3.add_edge(2, 0)
        print(self.w_3._map)
        self.w_3.add_edge(0, 3, 1)
        self.w_3.add_edge(1, 3, 2, 0)
        self.w_3.add_edge(2, 3, 0, 1)
        print(self.w_3._map)

    def test_get_next_edge(self):
        """Test get_next_edge()"""
        assert_equal((3,2), self.w_3.get_next_edge(3,1))
        assert_equal((0,1), self.w_3.get_next_edge(0,2))
        assert_equal((0,2), self.w_3.get_next_edge(0,3))

    def test_get_prev_edge(self):
        """Test get_prev_edge()"""
        assert_equal((3,0), self.w_3.get_prev_edge(3,1))
        assert_equal((0,1), self.w_3.get_prev_edge(0,3))
        assert_equal((0,3), self.w_3.get_prev_edge(0,2))