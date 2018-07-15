import networkx as nx
from collections import namedtuple, defaultdict


__all__ = ["combinatorial_embedding_to_pos"]


def combinatorial_embedding_to_pos(embedding):
    """Assigns every node a (x, y) position based on the given embedding

    Parameters
    ----------
    embedding : nx.PlanarEmbedding
        A combinatorial embedding that maps each node to a list of nodes,
        which defines the order of the outgoing edges (in clockwise order)

    Returns
    -------
    pos : dict
        Maps each node to a tuple that defines the (x, y) position
    """
    if len(embedding.nodes()) < 4:
        # Position the node in any triangle
        default_positions = [(0, 0), (2, 0), (1, 1)]
        pos = {}
        for i, v in enumerate(embedding.nodes()):
            pos[v] = default_positions[i]
        return pos

    embedding, outer_face = triangulate_embedding(embedding)

    # The following dicts map a node to another node
    # If a node is not in the key set it means that the node is not yet in G_k
    # If a node maps to Nil then the corresponding subtree does not exist
    left_t_child = {}
    right_t_child = {}

    # The following dicts map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = get_canonical_ordering(embedding, outer_face)
    print("Canonical ordering: ", node_list)
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
                                                      delta_x, vk)
        wp = contour_neighbor_data.wp
        wp1 = contour_neighbor_data.wp1
        wq = contour_neighbor_data.wq
        wq1 = contour_neighbor_data.wq1
        adds_mult_tri = contour_neighbor_data.adds_mult_tri
        delta_x_wp_wq = contour_neighbor_data.delta_x_wp_wq

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        # Adjust offsets
        delta_x[vk] = (-y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq])//2
        y_coordinate[vk] = (y_coordinate[wp] + delta_x_wp_wq +
                            y_coordinate[wq]) // 2
        delta_x[wq] = delta_x_wp_wq - delta_x[vk]
        if adds_mult_tri:
            delta_x[wp1] -= delta_x[vk]

        # Install v_k:
        right_t_child[wp] = vk
        right_t_child[vk] = wq
        if adds_mult_tri:
            left_t_child[vk] = wp1
            right_t_child[wq1] = Nil
        else:
            left_t_child[vk] = Nil

    # 2. Phase Set absolute positions
    pos = dict()
    pos[v1] = (0, y_coordinate[v1])
    remaining_nodes = [v1]
    while remaining_nodes:
        parent_node = remaining_nodes.pop()

        # Calculate position for left child
        set_absolute_position(parent_node, left_t_child,
                              remaining_nodes, delta_x, y_coordinate, pos)
        # Calculate position for right child
        set_absolute_position(parent_node, right_t_child,
                              remaining_nodes, delta_x, y_coordinate, pos)
    return pos


def set_absolute_position(parent, tree, remaining_nodes, delta_x, y_coordinate, pos):
    child = tree[parent]
    parent_node_x = pos[parent][0]
    if child is not Nil:
        # Calculate pos of child
        child_x = parent_node_x + delta_x[child]
        pos[child] = (child_x, y_coordinate[child])
        # Remember to calculate pos of its children
        remaining_nodes.append(child)


def get_canonical_ordering(embedding, outer_face):
    """Returns a canonical ordering of the nodes

    Parameters
    ----------
    embedding : nx.PlanarEmbedding
        The embedding is already fully triangulated
    outer_face : list
        The nodes on the outer face of the graph

    Returns
    -------
    node_list : list
        All nodes in the canonical ordering
    """
    v1 = outer_face[0]
    v2 = outer_face[1]
    chords = defaultdict(int)  # Maps nodes to the number of their chords
    marked_nodes = set()
    ready_to_pick = set(outer_face)

    # Initialize outer_face_ccw_nbr (do not include v1 -> v2)
    outer_face_ccw_nbr = {}
    prev_nbr = v2
    for idx in range(2, len(outer_face)):
        outer_face_ccw_nbr[prev_nbr] = outer_face[idx]
        prev_nbr = outer_face[idx]
    outer_face_ccw_nbr[prev_nbr] = v1

    # Initialize outer_face_cw_nbr (do not include v2 -> v1)
    outer_face_cw_nbr = {}
    prev_nbr = v1
    for idx in range(len(outer_face)-1, 0, -1):
        outer_face_cw_nbr[prev_nbr] = outer_face[idx]
        prev_nbr = outer_face[idx]

    def is_outer_face_nbr(v, w):
        if v not in outer_face_ccw_nbr:
            return outer_face_cw_nbr[v] == w
        if v not in outer_face_cw_nbr:
            return outer_face_ccw_nbr[v] == w
        return outer_face_ccw_nbr[v] == w or outer_face_cw_nbr[v] == w

    def is_on_outer_face(v):
        return v not in marked_nodes and (v in outer_face_ccw_nbr.keys() or
                                          v == v1 or v == v2)  # TODO: We can remove the v1 or the v2 comparison check which

    # Initialize number of chords
    for v in outer_face:
        for nbr in embedding.get_neighbors(v):
            if is_on_outer_face(nbr) and not is_outer_face_nbr(v, nbr):
                chords[v] += 1
                ready_to_pick.discard(v)

    # Initialize canonical_ordering
    canonical_ordering = [None]*len(embedding.nodes())
    canonical_ordering[0] = v1
    canonical_ordering[1] = v2
    ready_to_pick.discard(v1)
    ready_to_pick.discard(v2)

    for k in range(len(embedding.nodes())-1, 1, -1):
        # 1. Pick v from ready_to_pick
        v = ready_to_pick.pop()
        marked_nodes.add(v)
        canonical_ordering[k] = v

        # v has exactly two neighbors on the outer face (wp and wq)
        wp = None
        wq = None
        for nbr in embedding.get_neighbors(v):  # TODO: Break if wp and wq are found
            if nbr in marked_nodes:
                # Only consider nodes that are not yet removed
                continue
            if is_on_outer_face(nbr):
                # nbr is either wp or wq
                if nbr == v1:
                    wp = v1
                elif nbr == v2:
                    wq = v2
                else:
                    if outer_face_cw_nbr[nbr] == v:
                        # nbr is wp
                        wp = nbr
                    else:
                        # nbr is wq
                        wq = nbr

        # Obtain new nodes on outer face (neighbors of v from wp to wq)
        wp_wq = [wp]  # TODO: Save the new_outerface
        nbr = wp
        while nbr != wq:
            # Get next next neighbor (lies clockwise on the outer face)
            next_nbr = embedding.ccw_nbr[v][nbr]
            wp_wq.append(next_nbr)
            # Update outer face
            outer_face_cw_nbr[nbr] = next_nbr
            outer_face_ccw_nbr[next_nbr] = nbr
            # Move to next neighbor of v
            nbr = next_nbr

        if len(wp_wq) == 2:
            # There was a chord between wp and wq, decrease number of chords
            chords[wp] -= 1
            if chords[wp] == 0:
                ready_to_pick.add(wp)
            chords[wq] -= 1
            if chords[wq] == 0:
                ready_to_pick.add(wq)
        else:
            new_face_nodes = set(wp_wq[1:-1])
            for w in new_face_nodes:
                ready_to_pick.add(w)
                for nbr in embedding.get_neighbors(w):
                    if is_on_outer_face(nbr) and not is_outer_face_nbr(w, nbr):
                        chords[w] += 1
                        ready_to_pick.discard(w)
                        if nbr not in new_face_nodes:
                            # Also increase chord for the neighbor
                            chords[nbr] += 1
                            ready_to_pick.discard(nbr)

    return canonical_ordering


def get_contour_neighbors(right_t_child, embedding, delta_x, vk):
    """Returns sorted neighbors of v_k that are in C_k-1 (w_p, ..., w_q)

    # TODO: Reformulate this explanation
    Consider the graph G_(k-1), which is the subgraph induced by node_list[0:k].
    We can obtain the contour of it using the embedding and numerate the nodes
    according to their absolute x position: w_1, ..., w_m.
    We return all neighbors of the node v_k=node_list[k] that are in this contour
    and keep the order described above. We also return the indices p and q, such
    that w_p is the first neighbor in the contour and w_q the last neighbor.

    Note that because the graph is fully triangulated there are no nodes between
    w_p and w_q that are not also a neighbor of v_k. So we only need to consider
    the neighbors of v_k. We can find out if a neighbor of v_k lies on C_(n-1) by
    checking if an entry in right_t_child already exists.
    This way we can filter the neighbors of v_k such that only the ones in C_(n-1)
    remain. We then need to find the first and last node in x direction to get w_p
    and w_q. For i in [p, ..., q-1] it holds that right_t_child[w_i] is a neighbor
    of v_k, but for w_q this does not hold. So we have a way to determine w_q and
    because we know in which sequence the neighbors occur we directly have w_(q-1)
    w_p and w_(p+1).

    Parameters
    ----------
    right_t_child : dict
    embedding : nx.PlanarEmbedding
    delta_x : dict
    vk : node
    """
    # Calculate neighbors of v_k on the contour of G_(k-1)
    all_neighbors = embedding.get_neighbors(vk)
    contour_neighbors = list(filter(lambda w: w in right_t_child, all_neighbors))
    contour_neighbors_set = set(contour_neighbors)

    # Determine idx of w_q in contour_neighbors
    for q_idx in range(len(contour_neighbors)):
        if right_t_child[contour_neighbors[q_idx]] not in contour_neighbors_set:
            # contour_neighbors[idx] is w_q
            break

    # Determine all relevant contour neighbors
    wq = contour_neighbors[q_idx]
    wq1 = contour_neighbors[(q_idx + 1) % len(contour_neighbors)]
    wp = contour_neighbors[q_idx-1]
    wp1 = contour_neighbors[q_idx-2]  # Note that len(contour_neighbors) >= 2

    # Determine if v_k only adds multiple triangles
    adds_mult_tri = len(contour_neighbors) > 2

    # Calculate delta x(w_p, w_q)
    delta_x_wp_wq = 2  # +2 because the gaps are later stretched
    for idx in range(len(contour_neighbors)):
        if idx != (q_idx + 1) % len(contour_neighbors):  # Exclude w_p
            delta_x_wp_wq += delta_x[contour_neighbors[idx]]

    return ContourNeighborData(wp, wp1, wq1, wq, delta_x_wp_wq, adds_mult_tri)


def triangulate_face(embedding, v1, v2):
    """Triangulates the face given by half edge (v, w)"""
    _, v3 = embedding.next_face_half_edge(v1, v2)
    _, v4 = embedding.next_face_half_edge(v2, v3)
    if v1 == v2 or v1 == v3:
        # The component has less than 3 nodes
        return
    while v1 != v4:
        # Add edge if not already present on other side
        if embedding.has_edge(v1, v3):
            # Cannot triangulate at this position
            v1, v2, v3 = v2, v3, v4
        else:
            # Add edge for triangulation
            print("Add edge: ", v1, " ", v3)
            embedding.add_half_edge_cw(v1, v3, v2)
            embedding.add_half_edge_ccw(v3, v1, v2)
            v1, v2, v3 = v1, v3, v4
        # Get next node
        _, v4 = embedding.next_face_half_edge(v2, v3)


def triangulate_embedding(embedding):
    """Triangulates the embedding.

    Chooses the face that has the most vertices as the outer face. All other
    faces are triangulated by adding edges in a zig-zag way to keep the maximum
    degree low.

    The triangulated graph is 2-connected. This means on a cycle on the outer face
    a node cannot appear twice.

    Parameters
    ----------
    embedding : nx.PlanarEmbedding
        Maps each node to a list of neighbors in clockwise orientation. The input
        graph contains at least 3 nodes.

    Returns
    -------
    embedding : nx.PlanarEmbedding
        A new embedding containing all edges from the input embedding and the
        additional edges to internally triangulate the graph.

    start_triangle : tuple
        A tuple of 3 nodes (v1, v2, v3) that define a triangle in the graph.
        The edge (v1,v2) must lie on the outer face. When viewed with the edge
        (v1, v2) at the bottom, the node v1 lies to the left of v2.

    """
    if len(embedding.nodes()) <= 1:
        return embedding
    embedding = embedding.copy()
    G = embedding.get_graph()
    # Get a list with a node for each connected component
    component_nodes = [next(iter(x)) for x in nx.connected_components(G)]

    # 1. Make graph a single component (add edge between components)
    for i in range(len(component_nodes)-1):
        v1 = component_nodes[i]
        v2 = component_nodes[i+1]
        embedding.add_edge(v1, v2)

    # 2. Calculate faces, ensure 2-connectedness and determine outer face
    outer_face = []  # A face with the most number of nodes
    face_list = []
    edges_visited = set()
    for v in embedding.nodes():
        for w in embedding.get_neighbors(v):
            new_face = make_bi_connected(embedding, v, w, edges_visited)
            if new_face:
                # Found a new face
                print("Found face: ", new_face)
                face_list.append(new_face)
                if len(new_face) > len(outer_face):
                    outer_face = new_face

    print("Outer face ", outer_face)
    print("Face list", face_list)
    # 3. Triangulate internal faces
    for face in face_list:
        if face is outer_face:  # TODO: Test the output if th
            continue

        print("Current face: ", face)
        triangulate_face(embedding, face[0], face[1])

    return embedding, outer_face


def make_bi_connected(embedding, starting_node, outgoing_node, edges_counted):
    """Makes the face given by (starting_node, outgoing_node) 2-connected

    Parameters
    ----------
    embedding: nx.PlanarEmbedding
        The embedding that defines the faces
    starting_node : node
        A node on the face
    outgoing_node : node
        A node such that the half edge (starting_node, outgoing_node) belongs
        to the face
    edges_counted: set
        Set of all half-edges that belong to a face that has been counted
    Returns
    -------
    face_nodes: list
        A list of all nodes at the border of this face
    """

    # Check if the face has already been calculated
    if (starting_node, outgoing_node) in edges_counted:
        # This face was already counted
        return []
    edges_counted.add((starting_node, outgoing_node))

    # Add all edges to edges_counted which have this face to their left
    v1 = starting_node
    v2 = outgoing_node
    face_list = [starting_node]
    face_set = set(face_list)
    _, v3 = embedding.next_face_half_edge(v1, v2)
    while v2 != starting_node or v3 != outgoing_node:
        # cycle is not completed yet
        if v2 in face_set:
            print("Added biconnect edge: ", v1, " ", v3)
            embedding.add_half_edge_cw(v1, v3, v2)
            embedding.add_half_edge_ccw(v3, v1, v2)
            edges_counted.add((v2, v3))
            edges_counted.add((v3, v1))
            v2 = v1

        else:
            face_set.add(v2)
            face_list.append(v2)
        
        # set next edge
        v1 = v2
        v2, v3 = embedding.next_face_half_edge(v2, v3)

        # remember that this edge has been counted
        edges_counted.add((v1, v2))

    return face_list


ContourNeighborData = namedtuple('ContourNeighborData',
                                 ['wp', 'wp1', 'wq1', 'wq', 'delta_x_wp_wq',
                                  'adds_mult_tri'])


class Nil:
    """A class to represent that a node is not present

    We cannot use None, because None might be a node in the graph.
    """
    pass


def main():
    import matplotlib.pyplot as plt
    while True:
        embedding_data = {0: [1, 8], 1: [0, 3, 7], 2: [4, 7], 3: [1, 8], 4: [8, 2, 7, 5, 9], 5: [4], 6: [7], 7: [2, 1, 4, 6], 8: [3, 0, 4], 9: [4]}
        embedding = nx.PlanarEmbedding()
        embedding.set_data(embedding_data)
        G = embedding.get_graph()
        #embedding = None
        if not embedding:
            n = 10
            p = 0.9
            is_planar = False
            while not is_planar:
                G = nx.fast_gnp_random_graph(n, p)
                is_planar, embedding = nx.check_planarity(G)
                p /= 2
        print("Embedding ", embedding.get_data())
        pos = combinatorial_embedding_to_pos(embedding)
        print("Pos ", pos)
        nx.draw(G, pos)  # networkx draw()
        plt.draw()  # pyplot draw()
        plt.show()


if __name__ == '__main__':
    main()
