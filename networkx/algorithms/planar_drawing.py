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
    if len(embedding) < 4:
        pass
        # TODO: Special case for n < 4

    # The following dics map a node to another node
    left_t_child = {}
    right_t_child = {}

    # The following dics map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = get_canonical_ordering(embedding)

    # 1. Phase

    # Initialization
    v1, v2, v3 = node_list[0], node_list[1], node_list[2]

    delta_x[v1] = 0
    y_coordinate[v1] = 0
    right_t_child[v1] = v3
    left_t_child[v1] = Nil

    delta_x[v2] = 1
    y_coordinate[v2] = 0
    right_t_child[v2] = Nil
    left_t_child[v2] = Nil

    delta_x[v3] = 1
    y_coordinate[v3] = 1
    right_t_child[v3] = v2
    left_t_child[v3] = Nil

    for k in range(3, len(node_list)):
        neighbours, (p, q) = get_canonical_neighbors(node_list, k, embedding)
        wp, wp1, wq1, wq = neighbours[0], neighbours[1], neighbours[-2], neighbours[-1]
        vk = node_list[k]

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        # Adjust offsets
        delta_x_wp_wq = sum((delta_x[neighbours[i]] for i in
                             range(1, len(neighbours))))
        delta_x[vk] = (-y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq]) // 2
        y_coordinate[vk] = (y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq]) // 2
        delta_x[wq] = delta_x_wp_wq - delta_x[vk]
        if p + 1 != q:
            delta_x[wp1] -= delta_x[vk]

        # Install v_k:
        right_t_child[wp] = vk
        right_t_child[vk] = wq
        if p + 1 != q:
            left_t_child[wp] = vk
            right_t_child[wq1] = Nil
        else:
            left_t_child[vk] = Nil

    # 2. Phase
    accumulate_offsets(v1, 0, left_t_child, right_t_child, delta_x)


def get_canonical_ordering(embedding):
    """Returns a canonical ordering of the nodes

    # TODO: Maybe it might be of advantage to have not only the embedding but
    # the whole graph in this function (because of the dict structure).

    Parameters
    ----------
    embedding : dict

    Returns
    -------
    node_list : list
        All nodes in the canonical ordering
    """
    return [...]  # TODO: Implement (should return list of nodes)


def get_canonical_neighbors(node_list, k, embedding):
    """Returns sorted neighbors of v_k that are in C_k-1 (w_p, ..., w_q)
    TODO: Write docstring

    TODO: Is it more efficient to compute this once for all k?
    """
    return [...], (p, q)  # TODO: Implement (should return a list of nodes and a tuple (p, q))


def accumulate_offsets(vertex, delta, left_t_child, right_t_child, delta_x):
    """
    TODO: Write docstring
    """
    if vertex is Nil:
        return

    delta_x[vertex] += delta
    accumulate_offsets(left_t_child[vertex], delta_x[vertex],
                       left_t_child, right_t_child, delta_x)
    accumulate_offsets(right_t_child[vertex], delta_x[vertex],
                       left_t_child, right_t_child, delta_x)


class Nil:
    """A class to represent that a node is not present

    We cannot use None, because None might be a node in the graph.
    """
    pass
