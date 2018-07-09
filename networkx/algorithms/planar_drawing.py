import networkx as nx
from collections import namedtuple, defaultdict


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
        # Position the node in any triangle
        default_positions = [(0, 0), (2, 0), (1, 1)]
        pos = {}
        for i, v in enumerate(embedding):
            pos[v] = default_positions[i]
        return pos

    embedding, start_triangle = triangulate_embedding(embedding)

    # The following dicts map a node to another node
    # If a node is not in the key set it means that the node is not yet in G_k
    # If a node maps to Nil then the corresponding subtree does not exist
    left_t_child = {}
    right_t_child = {}

    # The following dicts map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = get_canonical_ordering(embedding, start_triangle)

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


def get_canonical_ordering(embedding, start_triangle):
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

    # Create dict holding for each node a dict mapping
    # neighbor vertices to index in the embedding list
    emb_idx = {}
    for u, l in embedding.items():
        emb_idx[u] = {}
        for idx, v in enumerate(l):
            emb_idx[u][v] = idx

    # Set v1, v2 and v3
    v1, v2, v3 = start_triangle

    # For not fully triangulated embeddings we might need to flip v1 and v2
    if embedding[v3][(emb_idx[v3][v2] + 1) % len(embedding[v3])] is not v1:
        # v1, v2, v3 is not a triangle
        v1, v2 = v2, v1
        v3 = embedding[v2][(emb_idx[v2][v1] + 1) % len(embedding[v2])]

    # Maintain a list for the result and a set for fast queries
    node_list = [v1, v2]
    node_set = set(node_list)

    # Remaining node stack
    insertable_nodes = [v3]

    # Obtain remaining order
    while len(node_list) != len(embedding):
        vk = insertable_nodes.pop()
        if vk not in node_set:
            # vk is the next node in the canonical ordering
            node_set.add(vk)
            node_list.append(vk)

            for v_next in embedding[vk]:
                idx = emb_idx[v_next][vk]
                if (embedding[v_next][(idx + 1) % len(embedding[v_next])] in node_set or
                        embedding[v_next][idx - 1] in node_set):
                    insertable_nodes.append(v_next)

    return node_list


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
    """
    # Calculate neighbors of v_k on the contour of G_(k-1)
    all_neighbors = embedding[vk]
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



def get_face(embedding, edge, marked_edges):
    if edge in marked_edges: #face was already found once
        return None
    v,u = edge

    face = []
    face.append(v)
    first_node = v
    edge_start=v
    edge_end = u

    while edge_end != v:
        face.append(edge_end) #TODO: check if it already contains edge_end then additional measures for 2 connectedness are necessary
        nghbrs = embedding[edge_end]
        for idx,n in enumerate(nghbrs):
            if n == edge_start:
                next_node = nghbrs[(idx+1) % len(nghbrs)]
                break
        edge_start=edge_end
        edge_end=next_node
        marked_edges.add((edge_start,edge_end))

    return face


def triangulate_embedding(embedding):
    """Triangulates the embedding.

    Chooses the face that has the most vertices as the outer face. All other
    faces are triangulated by adding edges in a zig-zag way to keep the maximum
    degree low.

    The triangulated graph is 2-connected. This means on a cycle on the outer face
    a node cannot appear twice.

    Parameters
    ----------
    embedding : dict
        Maps each node to a list of neighbors in clockwise orientation. The input
        graph contains at least 3 nodes.

    Returns
    -------
    embedding : dict
        A new embedding containing all edges from the input embedding and the
        additional edges to internally triangulate the graph.

    start_triangle : tuple
        A tuple of 3 nodes (v1, v2, v3) that define a triangle in the graph.
        The edge (v1,v2) must lie on the outer face. When viewed with the edge
        (v1, v2) at the bottom, the node v1 lies to the left of v2.

    """
     # TODO: Finish Implementation
    if len(embedding) <= 1:
        return embedding

    # 1. Transform embedding datastructure (and obtain components)
    left_neighbor = defaultdict(dict)  # Maps a node v to a dict that maps a neighbor of v to the neighbor left of it
    right_neighbor = defaultdict(dict)  # Similar to the above but for the right neighbor
    # TODO: Component tracking

    #Get components
    component_nodes = []  # A node for each component in the graph

    #Just DFS to get the components
    visited_nodes = set()
    #Think about merging this loop into the one below

    for v, neighbor_list in embedding.items():
        if v not in visited_nodes:
           component_nodes.append(v)
           visited_nodes.add(v)
           to_visit=[x for x in neighbor_list]
           while len(to_visit) > 0:
               n=to_visit.pop()
               if n not in visited_nodes:
                   visited_nodes.add(n)
                   to_visit.extend(embedding[n])

    for v, neighbor_list in embedding.items():
        for i in range(len(neighbor_list)):
            left_neighbor1 = neighbor_list[i-1]
            right_neighbor1 = neighbor_list[i]
            left_neighbor[v][right_neighbor1] = left_neighbor1
            right_neighbor[v][left_neighbor1] = right_neighbor1

    # 2. Make graph a single component (add edge between components)
    for i in range(len(component_nodes)-1):
        v1 = component_nodes[i]
        v2 = component_nodes[i+1]
        add_half_edge(v1, v2, left_neighbor, right_neighbor)
        add_half_edge(v2, v1, left_neighbor, right_neighbor)

    # 3. Calculate faces, ensure 2-connected and determine outer face
    faces = []
    marked_edges = set()


    #edges=[]
    #     for v,neighbor_list in embedding.items():
    #    edges.extend(((v,x) for x in neighbor_list))

    for v,neighbor_list in embedding.items():
        for u in neighbor_list:
            face = get_face(embedding,(v,u),marked_edges)
            if face is not None:
                faces.append(face)

    outer_face =  max(faces, key=lambda x:len(x))  # A face with the most number of nodes

    # 4. Triangulate internal faces (in a sig-zack fashion)

    for face in faces:
        if face != outer_face:
            edge_start=0
            edge_end=0;
            dir=1
            for i in range(1,len(face)):
                edge_end = (edge_end + i*dir) % len(face) #Zig zag pattern
                add_half_edge(edge_start, edge_end, left_neighbor, right_neighbor)
                add_half_edge(v1, v2, edge_start, edge_end)
                edge_start = edge_end



    # 5. Transform embedding datastructure back
    new_embedding = dict()
    for v in embedding:
        start_node = next(iter(right_neighbor[v]))
        current_node = start_node
        neighbor_list = list()
        while len(neighbor_list) == 0 or start_node != current_node:
            neighbor_list.append(current_node)
            current_node = right_neighbor[v][current_node]
        new_embedding[v] = neighbor_list

    # 6. Choose start_triangle
    v1 = outer_face[0]
    v2 = outer_face[1]
    v3 = left_neighbor[v1][v2]
    return new_embedding, (v1, v2, v3)

ContourNeighborData = namedtuple('ContourNeighborData',
                                 ['wp', 'wp1', 'wq1', 'wq', 'delta_x_wp_wq',
                                  'adds_mult_tri'])


class Nil:
    """A class to represent that a node is not present

    We cannot use None, because None might be a node in the graph.
    """
    pass
