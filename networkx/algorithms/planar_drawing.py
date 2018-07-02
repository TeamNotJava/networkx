import networkx as nx
from collections import namedtuple


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
        raise NotImplementedError
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
        vk = node_list[k]
        contour_neighbor_data = get_contour_neighbors(right_t_child, embedding,
                                                      delta_x, v1, vk)
        wp = contour_neighbor_data.wp
        wp1 = contour_neighbor_data.wp1
        wq = contour_neighbor_data.wq
        wq1 = contour_neighbor_data.wq1
        p = contour_neighbor_data.p
        q = contour_neighbor_data.q
        delta_x_wp_wq = contour_neighbor_data.delta_x_wp_wq

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        # Adjust offsets
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
            left_t_child[vk] = wp1
            right_t_child[wq1] = Nil
        else:
            left_t_child[vk] = Nil

    # 2. Phase
    accumulate_offsets(v1, 0, left_t_child, right_t_child, delta_x)

    # 3. Phase: Calculate absolute positions
    pos = dict()
    pos[v1] = (0, y_coordinate[v1])
    remaining_nodes = [v1]
    while remaining_nodes:
        parent_node = remaining_nodes.pop()
        parent_node_x = pos[parent_node][0]

        left_child = left_t_child[parent_node]
        if left_child is not Nil:
            # Calculate pos of left child
            left_child_x = parent_node_x + delta_x[left_child]
            pos[left_child] = (left_child_x, y_coordinate[left_child])
            # Remember to calculate pos of its children
            remaining_nodes.append(left_child)

        right_child = right_t_child[parent_node]
        if right_child is not Nil:
            # Calculate pos of right child
            right_child_x = parent_node_x + delta_x[right_child]
            pos[right_child] = (right_child_x, y_coordinate[right_child])
            # Remember to calculate pos of its children
            remaining_nodes.append(right_child)


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


def get_contour_neighbors(right_t_child, embedding, delta_x, v1, vk):
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
    TODO: Not sure if w_(p+1) and w_(q-1) must also be neighbors of v_k
    TODO: Not sure if delta_x_wp_wq should only be computed over neighbors of v_k
    """
    contour_node = v1
    neighbor_set = set(embedding[vk])
    p, q = None, None
    wp, wp1, wq1, wp = None, None, None, None
    delta_x_wp_wq = 2  # +2 because the gaps are later stretched
    delta_x_wp_wq_temp = 0
    idx = 0

    while True:
        if p is not None:  # idx > p
            delta_x_wp_wq_temp += delta_x[contour_node]

        if contour_node in neighbor_set:
            if wp is None:
                # The first contour_node that is a neighbor of vk
                p = idx
                wp = contour_node
            # It might be the last contour_node that is a neighbor of vk
            q = idx
            wq = contour_node
            wq1 = maybe_wq1  # Set to previously encountered contour_node

            delta_x_wp_wq += delta_x_wp_wq_temp
            delta_x_wp_wq_temp = 0

        maybe_wq1 = contour_node
        if idx + 1 == p:
            wp1 = contour_node
        # Get the next contour_node:
        idx += 1
        if right_t_child[contour_node] is not Nil:
            contour_node = right_t_child[contour_node]
        else:
            break

    return ContourNeighborData(wp, wp1, wq1, wq, delta_x_wp_wq, p, q)


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


ContourNeighborData = namedtuple('ContourNeighborData',
                                 ['wp', 'wp1', 'wq1', 'wp', 'delta_x_wp_wq',
                                  'p', 'q'])


class Nil:
    """A class to represent that a node is not present

    We cannot use None, because None might be a node in the graph.
    """
    pass
