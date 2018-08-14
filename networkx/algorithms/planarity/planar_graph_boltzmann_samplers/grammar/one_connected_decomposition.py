# -*- coding: utf-8 -*-
#    Copyright (C) 2018 by
#    Marta Grobelna <marta.grobelna@rwth-aachen.de>
#    Petre Petrov <petrepp4@gmail.com>
#    Rudi Floren <rudi.floren@gmail.com>
#    Tobias Winkler <tobias.winkler1@rwth-aachen.de>
#    All rights reserved.
#    BSD license.
#
# Authors:  Marta Grobelna <marta.grobelna@rwth-aachen.de>
#           Petre Petrov <petrepp4@gmail.com>
#           Rudi Floren <rudi.floren@gmail.com>
#           Tobias Winkler <tobias.winkler1@rwth-aachen.de>

from ...boltzmann_framework.evaluation_oracle import EvaluationOracle
from ...boltzmann_framework.generic_samplers import *
from ...boltzmann_framework.decomposition_grammar import AliasSampler, DecompositionGrammar
from ...boltzmann_framework.generic_samplers import BoltzmannSamplerBase

from ..combinatorial_classes.one_connected_graph import OneConnectedPlanarGraph
from .grammar_utils import underive
from .two_connected_decomposition import two_connected_graph_grammar


class Merger(DefaultBuilder):
    """
    Merges a set of l-derived graphs at their marked vertices.
    """
    def set(self, graphs):
        # Merge a set of l-derived one-connected planar graphs at their marked vertices.
        # If the set is empty, return a single-node graph.
        if len(graphs) is 0:
            g = OneConnectedPlanarGraph()
            return LDerivedClass(OneConnectedPlanarGraph(), g.half_edge)
        result = graphs.pop()
        for g in graphs:
            result.marked_atom.insert_all(g.marked_atom)
        assert isinstance(result, LDerivedClass)
        return result


def merge(prod):
    """Merges l-derived one-connected graphs at their marked vertices"""
    # lhs is a bi-derived connected and rhs a derived connected.
    lhs = prod.first
    rhs = prod.second
    if not lhs.base_class_object.marked_atom.is_trivial:
        rhs.marked_atom.insert_all(lhs.base_class_object.marked_atom)
    return lhs


def subs_marked_vertex(decomp):
    # decomp is of form (G_1_dx + L * G_1_dx_dx) * G_2_dx_dx.
    if isinstance(decomp.first, LDerivedClass):
        plug_in_he = decomp.first.marked_atom
    else:
        plug_in_he = decomp.first.second.base_class_object.marked_atom # or ... .base_class_object.marked_atom ?
    if not plug_in_he.is_trivial:
        decomp.second.marked_atom.insert_all(plug_in_he)
    return decomp.second


def rej_to_G_1(g):
    return bern(1 / (g.l_size + 1))


def one_connected_graph_grammar():
    """Constructs the grammar for connected planar graphs.

    Returns
    -------
    DecompositionGrammar
        The grammar for sampling from G_1_dx and G_1_dx_dx.
    """

    # Some shortcuts to make the grammar more readable.
    L = LAtomSampler
    G_2_dx = AliasSampler('G_2_dx')
    G_2_dx_dx = AliasSampler('G_2_dx_dx')
    G_1_dx = AliasSampler('G_1_dx')
    G_1_dx_dx = AliasSampler('G_1_dx_dx')
    Set = SetSampler
    LSubs = LSubsSampler
    Bij = BijectionSampler
    Rej = RejectionSampler

    grammar = DecompositionGrammar()
    grammar.rules = two_connected_graph_grammar().rules
    grammar.rules = {

        # 1 connected planar graphs

        'G_1_dx': Set(0, LSubs(G_2_dx, L() * G_1_dx)),

        'G_1_dx_dx':
            Bij(
                Bij((G_1_dx + L() * G_1_dx_dx) * LSubs(G_2_dx_dx, L() * G_1_dx), subs_marked_vertex) * G_1_dx,
            merge),

        'G_1': Bij(Rej(G_1_dx, rej_to_G_1), underive)  # See lemma 15.

    }
    grammar.set_builder(['G_1_dx'], Merger())

    return grammar
