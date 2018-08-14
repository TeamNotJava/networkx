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

from ...boltzmann_framework.generic_samplers import *
from ...boltzmann_framework.decomposition_grammar import DecompositionGrammar, AliasSampler
from ...boltzmann_framework.evaluation_oracle import EvaluationOracle
from ...boltzmann_framework.generic_samplers import BoltzmannSamplerBase

from ..bijections.closure import Closure
from ..combinatorial_classes import BinaryTree
from .binary_tree_decomposition import binary_tree_grammar


def closure(binary_tree):
    """To be used as bijection in the grammar."""
    if isinstance(binary_tree, LDerivedClass):
        binary_tree = binary_tree.base_class_object
    assert isinstance(binary_tree, BinaryTree), binary_tree
    dissection = Closure().closure(binary_tree)
    return dissection


def add_random_root_edge(decomp):
    """From ((L, U), dissection) or (U, dissection) to IrreducibleDissection."""
    dissection = decomp.second
    dissection.root_at_random_hexagonal_edge()
    return dissection


def is_admissible(dissection):
    """Admissibility check for usage in the grammar."""
    return dissection.is_admissible


def irreducible_dissection_grammar():
    """Builds the dissection grammar. Must still be initialized with init().

    Returns
    -------
    DecompositionGrammar
        The grammar for sampling from J_a and J_a_dx.
    """

    # Some shorthands to keep the grammar readable.
    L = LAtomSampler
    U = UAtomSampler
    K = AliasSampler('K')
    K_dx = AliasSampler('K_dx')
    I = AliasSampler('I')
    I_dx = AliasSampler('I_dx')
    J = AliasSampler('J')
    J_dx = AliasSampler('J_dx')
    Bij = BijectionSampler
    Rej = RejectionSampler

    grammar = DecompositionGrammar()
    # This grammar depends on the binary tree grammar so we add it.
    grammar.rules = binary_tree_grammar().rules

    grammar.rules = {

        'I': Bij(K, closure),

        'I_dx': Bij(K_dx, closure),

        'J': Bij(3*L()*U()*I, add_random_root_edge),

        'J_dx': Bij(3*U()*I + 3*L()*U()*I_dx, add_random_root_edge),

        'J_a': Rej(J, is_admissible),

        'J_a_dx': Rej(J_dx, is_admissible),

    }
    return grammar

