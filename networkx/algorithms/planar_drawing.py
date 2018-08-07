import networkx as nx
from collections import namedtuple, defaultdict


__all__ = ["combinatorial_embedding_to_pos"]


def combinatorial_embedding_to_pos(embedding, fully_triangulate=False):
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

    embedding, outer_face = triangulate_embedding(embedding, fully_triangulate)

    # The following dicts map a node to another node
    # If a node is not in the key set it means that the node is not yet in G_k
    # If a node maps to None then the corresponding subtree does not exist
    left_t_child = {}
    right_t_child = {}

    # The following dicts map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = get_canonical_ordering(embedding, outer_face)

    # 1. Phase

    # Initialization
    v1, v2, v3 = node_list[0][0], node_list[1][0], node_list[2][0]

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
        vk, contour_neighbors = node_list[k]
        wp = contour_neighbors[0]
        wp1 = contour_neighbors[1]
        wq = contour_neighbors[-1]
        wq1 = contour_neighbors[-2]
        adds_mult_tri = len(contour_neighbors) > 2

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        delta_x_wp_wq = sum((delta_x[x] for x in contour_neighbors[1:]))

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
            right_t_child[wq1] = None
        else:
            left_t_child[vk] = None

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
    if child is not None:
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
        for nbr in embedding.neighbors_cw_order(v):
            if is_on_outer_face(nbr) and not is_outer_face_nbr(v, nbr):
                chords[v] += 1
                ready_to_pick.discard(v)

    # Initialize canonical_ordering
    canonical_ordering = [None]*len(embedding.nodes())  # type: object
    canonical_ordering[0] = (v1, [])
    canonical_ordering[1] = (v2, [])
    ready_to_pick.discard(v1)
    ready_to_pick.discard(v2)

    for k in range(len(embedding.nodes())-1, 1, -1):
        # 1. Pick v from ready_to_pick
        v = ready_to_pick.pop()
        marked_nodes.add(v)

        # v has exactly two neighbors on the outer face (wp and wq)
        wp = None
        wq = None
        for nbr in embedding.neighbors_cw_order(v):  # TODO: Break if wp and wq are found
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
            if wp is not None and wq is not None:
                # We don't need to iterate any further
                break

        # Obtain new nodes on outer face (neighbors of v from wp to wq)
        wp_wq = [wp]
        nbr = wp
        while nbr != wq:
            # Get next next neighbor (lies clockwise on the outer face)
            next_nbr = embedding[v][nbr]['ccw']
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
            # Update all chords involving w_(p+1) to w_(q-1)
            new_face_nodes = set(wp_wq[1:-1])
            for w in new_face_nodes:
                # If we do not find a chord for w later we can pick it next
                ready_to_pick.add(w)
                for nbr in embedding.neighbors_cw_order(w):
                    if is_on_outer_face(nbr) and not is_outer_face_nbr(w, nbr):
                        # There is a chord involving w
                        chords[w] += 1
                        ready_to_pick.discard(w)
                        if nbr not in new_face_nodes:
                            # Also increase chord for the neighbor
                            # We only iterator over new_face_nodes
                            chords[nbr] += 1
                            ready_to_pick.discard(nbr)
        # Set the canonical ordering node and the list of contour neighbors
        canonical_ordering[k] = (v, wp_wq)

    return canonical_ordering


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
            embedding.add_half_edge_cw(v1, v3, v2)
            embedding.add_half_edge_ccw(v3, v1, v2)
            v1, v2, v3 = v1, v3, v4
        # Get next node
        _, v4 = embedding.next_face_half_edge(v2, v3)


def triangulate_embedding(embedding, fully_triangulate=False):
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
    embedding = nx.PlanarEmbedding(embedding)

    # Get a list with a node for each connected component
    component_nodes = [next(iter(x)) for x in nx.connected_components(embedding)]

    # 1. Make graph a single component (add edge between components)
    for i in range(len(component_nodes)-1):
        v1 = component_nodes[i]
        v2 = component_nodes[i+1]
        embedding.connect_components(v1, v2)

    # 2. Calculate faces, ensure 2-connectedness and determine outer face
    outer_face = []  # A face with the most number of nodes
    face_list = []
    edges_visited = set()
    for v in embedding.nodes():
        for w in embedding.neighbors_cw_order(v):
            new_face = make_bi_connected(embedding, v, w, edges_visited)
            if new_face:
                # Found a new face
                face_list.append(new_face)
                if len(new_face) > len(outer_face):
                    outer_face = new_face

    # 3. Triangulate internal faces
    for face in face_list:
        if face is outer_face and not fully_triangulate:
            # Don't triangulate this face
            continue
        triangulate_face(embedding, face[0], face[1])

    if fully_triangulate:
        v1 = outer_face[0]
        v2 = outer_face[1]
        v3 = embedding[v2][v1]['ccw']
        outer_face = [v1, v2, v3]

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



def main():
    import matplotlib.pyplot as plt

    while True:
        n = 50
        p = 1.0
        is_planar = False
        while not is_planar:
            G = nx.fast_gnp_random_graph(n, p)
            is_planar, embedding = nx.check_planarity(G)
            p *= 0.99
        print("Embedding: ", embedding.get_data())
        print("Displaying not fully triangulated drawing")
        plt.subplot(1, 2, 1)
        nx.draw_planar(embedding, node_size=2)


        pos = combinatorial_embedding_to_pos(embedding, fully_triangulate=True)
        print("Displaying fully triangulated drawing")
        plt.subplot(1, 2, 2)
        nx.draw(G, pos, node_size=2)

        plt.show()


if __name__ == '__main__':
    main()
