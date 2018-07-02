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

    embedding = triangulate_embedding(embedding)

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
        contour_neighbors, p, q = get_contour_neighbors(right_t_child,
                                                        embedding, v1, vk)
        wp, wp1 = contour_neighbors[0], contour_neighbors[1]
        wq1, wq = contour_neighbors[-2], contour_neighbors[-1]
        vk = node_list[k]

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        # Adjust offsets
        delta_x_wp_wq = sum((delta_x[contour_neighbors[i]] for i in
                             range(1, len(contour_neighbors))))
        delta_x[vk] = (-y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq])//2
        y_coordinate[vk] = (y_coordinate[wp] + delta_x_wp_wq +
                            y_coordinate[wq]) // 2
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

    # 3. Phase: Calculate absolute positions
    # TODO


def get_canonical_ordering(embedding):
    """Returns a canonical ordering of the nodes

    Parameters
    ----------
    embedding : dict
        The embedding is already fully triangulated

    Returns
    -------
    node_list : list
        All nodes in the canonical ordering
    """
    return [...]  # TODO: Implement (should return list of nodes)


def get_contour_neighbors(right_t_child, embedding, v1, vk):
    """Returns sorted neighbors of v_k that are in C_k-1 (w_p, ..., w_q)

    Consider the graph G_(k-1), which is the subgraph induced by node_list[0:k].
    We can obtain the contour of it using the embedding and numerate the nodes
    according to their absolute x position: w_1, ..., w_m.
    We return all neighbors of the node v_k=node_list[k] that are in this contour
    and keep the order described above. We also return the indices p and q, such
    that w_p is the first neighbor in the contour and w_q the last neighbor.

    Travers the tree T by the use of right_t_child only along the right to
    obtain the contour nodes.
    TODO: Check complexity of this implementation
    TODO: Is this implementation actually correct?
    """
    contour_node = v1
    neighbor_set = set(embedding[vk])
    contour_neighbors = []
    p, q = None, None
    idx = 0
    while True:
        if contour_node in neighbor_set:
            # The contour node is a neighbor of v_k
            contour_neighbors.append(contour_node)
            q = idx  # q is the last index, update it every time
            if p is None:
                p = idx  # p is the first index, only set on first occurrence
        idx += 1
        if contour_node in right_t_child and right_t_child[contour_node] is not Nil:
            contour_node = right_t_child[contour_node]
        else:
            break

    return contour_neighbors, p, q


def triangulate_embedding(embedding):
    """Triangulates the embedding.

    Adds edges to the embedding until all faces are triangles.
    """
    # TODO: Implement
    return new_embedding


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
