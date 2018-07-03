from .generic_classes import *


class CombinatorialClassBuilder:
    """
    Interface for objects that build combinatorial classes.
    """

    def zero_atom(self):
        raise NotImplementedError

    def l_atom(self):
        raise NotImplementedError

    def u_atom(self):
        raise NotImplementedError

    def product(self, lhs, rhs):
        raise NotImplementedError

    def set(self, elements):
        raise NotImplementedError


class DefaultCombinatorialClassBuilder(CombinatorialClassBuilder):
    """
    Builds the generic objects.
    """

    def zero_atom(self):
        return ZeroAtomClass()

    def l_atom(self):
        return LAtomClass()

    def u_atom(self):
        return UAtomClass()

    def product(self, lhs, rhs):
        return ProdClass(lhs, rhs)

    def set(self, elements):
        return SetClass(elements)

class BinaryTreeBuilder(CombinatorialClassBuilder):
    """
    Builds bicolored binary trees.
    """

    def zero_atom(self):
        pass

    def l_atom(self):
        pass

    def u_atom(self):
        pass

    def product(self, lhs, rhs):
        pass

    def set(self, elements):
        pass
    