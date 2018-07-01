import networkx as nx


__all__ = ["combinatorial_embedding_to_pos"]


def combinatorial_embedding_to_pos(embedding):
    """Assigns every node a (x, y) position based on the given embedding

    Parameters
    ----------
    embedding : dict
        A combinatorial embedding that maps each node to a list of nodes,
        which defines the order of the outgoing edges (in clockwise order)

    Returns
    -------
    pos : dict
        Maps each node to a tuple that defines the (x, y) position
    """
    # TODO: What to do if graph contains node that has value None?
    # TODO: Special case for n < 4

    # The following dics map a node to another node
    left_t_child = {}
    right_t_child = {}

    # The following dics map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = [...]  # TODO: Canonical Ordering

    # 1. Phase

    # Initialization
    v1, v2, v3 = node_list[0], node_list[1], node_list[2]

    delta_x[v1] = 0
    y_coordinate[v1] = 0
    right_t_child[v1] = v3
    left_t_child[v1] = None

    delta_x[v2] = 1
    y_coordinate[v2] = 0
    right_t_child[v2] = None
    left_t_child[v2] = None

    delta_x[v3] = 1
    y_coordinate[v3] = 1
    right_t_child[v3] = v2
    left_t_child[v3] = None

    for k in range(3, len(node_list)):
        pass  # TODO: Implement

    # TODO: Finish Phase 1

    # 2. Phase
    accumulate_offsets(v1, 0, left_t_child, right_t_child, delta_x)


def accumulate_offsets(vertex, delta, left_t_child, right_t_child, delta_x):
    """
    TODO: Write docstring
    """
    if vertex is None:
        return

    delta_x[vertex] += delta
    accumulate_offsets(left_t_child[vertex], delta_x[vertex],
                       left_t_child, right_t_child, delta_x)
    accumulate_offsets(right_t_child[vertex], delta_x[vertex],
                       left_t_child, right_t_child, delta_x)
