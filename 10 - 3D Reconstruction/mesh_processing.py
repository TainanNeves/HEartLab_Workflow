import numpy as np
from tqdm import tqdm
from scipy.sparse import csr_array, issparse


#Esta função calcula a distância dos vértices de uma malha 3D para uma região específica
def mesh_distance_from_region(faces, reg_idx):
    connectivity = get_connectivity(faces,)
    size = faces.max() +1
    distances = np.zeros((size)) -1
    distances[reg_idx] = 0
    d = 1
    distances_old = float('nan')
    while np.any(distances==-1):
        if np.all(distances == distances_old):
            # print("Something went wrong while calculating nodal distance")
            distances[distances == -1] = d+1
            break
        distances_old = distances.copy()
        distances[(distances == -1)*connectivity[distances==(d-1)].sum(axis=0).astype(bool)] = d
        d += 1

    return distances

def get_connectivity(triangles, autoindex=True):
    """Get 1st and 2nd order connectivity matrix from triangles."""
    size = triangles.max() + 1
    connectivity = np.zeros((size, size), bool)
    for i, j, k in triangles:
        connectivity[i,j] = True
        connectivity[i,k] = True
        connectivity[j,i] = True
        connectivity[k,i] = True
        connectivity[j,k] = True
        connectivity[k,j] = True
    
    if autoindex:
        connectivity[np.arange(size), np.arange(size)] = True
    else:
        connectivity[np.arange(size), np.arange(size)] = False
    
    return connectivity

def find_HR_idx(HR_vertices, LR_vertices):
    idx = np.zeros(LR_vertices.shape[0], int)
    for i in range(LR_vertices.shape[0]):
        dist = HR_vertices - LR_vertices[(i,),:]
        idx[i] = np.argmin(np.linalg.norm(dist, axis=-1))
        
    return idx

def mesh_distance(faces, max_distance=-1, start_distances=None):
    if max_distance <= 0:
        max_distance = float('inf')
    connectivity = get_connectivity(faces,)
    size = faces.max() +1
    if start_distances is None:
        distances = np.eye(size) -1
    else:
        distances = start_distances.copy()
        distances[distances == float('inf')] = -1
        distances[distances < 0] = -1

    conn_i = connectivity.copy()
    
    d = distances.max() + 1
    while d <= max_distance:
        if d > 1:
            for i in range(size):
                conn_i[i] = connectivity[conn_i[i],:].any(axis=0)
        
        if not np.any((distances == -1) & conn_i):
            break                
        
        distances[(distances == -1) & conn_i] = d
        d += 1

    distances[distances == -1] = float('inf')
    return distances

def mesh_distance2(faces, max_distance=-1, start_distances=None):
    if max_distance <= 0:
        max_distance = faces.max() + 1
    connectivity = get_connectivity(faces,)
    size = faces.max() +1
    conn = connectivity.copy()

    if start_distances is None:
        distances = conn * 2 - np.eye(size) - 1

    else:
        distances = start_distances.copy()

    distances[distances < 0] = max_distance + 1
        
    
    d = np.sort(np.unique(distances))[-2] + 1
    while d <= max_distance:
        if d < 5:
            if d > 1:
                for i in range(size):
                    conn[i] = connectivity[conn[i],:].any(axis=0)
            
            if not np.any((distances >= max_distance) & conn):
                break                
            
            distances[(distances >= max_distance) & conn] = d
            d += 1

        else:
            distances_i = distances.copy()
            for i in range(size):
                distances_i[i, distances[i] <= max_distance] = np.min(distances_i[i, distances[i] <= max_distance].reshape(-1, 1) + distances[distances[i] <= max_distance], axis=1)
            
            distances = distances_i
            distances[distances > max_distance] = max_distance + 1

            d = np.unique(distances[distances <= max_distance]).max() + 1

    distances[distances > max_distance] = float('inf')
    return distances

#%%
#Código para implementar uma operação de erosão em uma máscara 3D binária usando conectividade
def erodemorph3d(mask, connectivity, r=1):
    for _ in range(r):
        mask = ~(connectivity[~mask].sum(axis=0) >0)
    
    return mask

#Esta função implementa uma operação de dilatação em uma máscara 3D binária usando conectividade
def dilatemorph3d(mask, connectivity, r=1):
    for _ in range(r):
        mask = connectivity[mask].sum(axis=0) >0
    
    return mask






#%%
#Esta função é utilizada para "mascarar" (selecionar) um subconjunto de triângulos de alta resolução (HR) com base em uma máscara ou índices fornecidos.
def mask_triangles(mask_or_inds, HR_triangles):
    
    if mask_or_inds.size > 1 and mask_or_inds.max() <= 1:
        inds = np.where(mask_or_inds == 1)[0]
        
    else:
        inds = mask_or_inds
    
    isin_mask = np.zeros_like(HR_triangles, bool)
    
    for ind in inds:
        isin_mask |= HR_triangles == ind
        
    LR_triangles = HR_triangles[np.logical_and(isin_mask, axis=1)]
    
    # for i, ind in enumerate(np.unique(LR_triangles)):
    #     LR_triangles[LR_triangles == ind] = i
    
    return LR_triangles

#Esta função é utilizada para encontrar as arestas de fronteira de uma malha 3D representada pelos triângulos fornecidos.
def border_edges(triangles):
    triangles.sort(axis=1)
    
    edges = np.concatenate([triangles[:,(0, 1),], triangles[:,(0, 2),].tolist(), triangles[:,(1, 2),]], axis=0)
    
    border = np.zeros(edges.shape[0], bool)
    taked = np.zeros(edges.shape[0], bool)
    
    for i in range(edges.shape[0]):
        if taked[i]:
            continue
        
        taked[i] = True
        
        same_edge = np.sum(edges[i] == edges[~taked], axis=1) == 2
        
        if same_edge.sum() > 0:
            taked[~taked] += same_edge
            
        else:
            border[i] = True
            
    
    return edges[border]

def get_scalar_connectivity(triangles):
    """Get 1st and 2nd order connectivity matrix from triangles."""
    size = triangles.max() + 1
    connectivity = np.zeros((size, size))
    for i, j, k in triangles:
        connectivity[i,j] += 1 
        connectivity[i,k] += 1
        connectivity[j,i] += 1
        connectivity[k,i] += 1
        connectivity[j,k] += 1
        connectivity[k,j] += 1
    
    return connectivity

def generate_edge_center_list(conn1, conn2):
    triplets = []
    for i, j in zip(*np.where(conn1)):
        if j > i:
            for k in np.where(conn2[[i, j],].all(axis=0))[0]:
                triplets.append([i, j, k])

    return triplets

def get_edges_from_faces(faces):
    edges = np.vstack([faces[:,[0, 1]], faces[:,[0, 2]], faces[:,[1, 2]]])
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    return edges

def mesh_neighbor_spatial_filter(signals, iterations=1, lamb=0.8, **kwargs):
    connectivity = kwargs.get("connectivity", get_connectivity(kwargs.get("faces", [0]), autoindex=False))
    connectivity[np.diag(connectivity)] = False
    for iter in range(iterations):
        signals_i = signals.copy()
        for i in range(signals.shape[0]):
            signals_i[i] = signals[i] * (1 - lamb) + signals[connectivity[i]].mean(axis=0) * lamb
        
        signals = signals_i
    
    return signals


def mesh_laplacian(vertex, face):
    nvertex = vertex.shape[0]
    nface = face.shape[0]

    print(f'MESH_LAPLACIAN: Calc Laplacian matrix for {nvertex} vertices... ', end="")
    
    # the matrix 'edge' is the connectivity of all vertices
    edge = np.zeros((nvertex, nvertex))
    
    for i in range(nface):
        # compute the length of all triangle edges (Diff is [3x3])
        Diff = vertex[face[i, [0, 1, 2]], :] - vertex[face[i, [1, 2, 0]], :]
        Norm = np.sqrt(np.sum(Diff**2, axis=1))

        edge[face[i, 0], face[i, 1]] = Norm[0]
        edge[face[i, 1], face[i, 2]] = Norm[1]
        edge[face[i, 2], face[i, 0]] = Norm[2]

        # make sure that all edges are symmetric
        edge[face[i, 1], face[i, 0]] = Norm[0]
        edge[face[i, 2], face[i, 1]] = Norm[1]
        edge[face[i, 0], face[i, 2]] = Norm[2]

    # Using edge to identify nearest vertices, calculate
    # the Laplacian for an irregular mesh
    lap = np.zeros((nvertex, nvertex))
    
    for i in range(nvertex):
        k = np.where(edge[i, :])[0]  # the indices of the neighbours
        ni = len(k)  # the number of neighbours
        hi = np.mean(edge[i, k])  # the average distance to the neighbours
        invhi = np.mean(1. / edge[i, k])  # the average inverse distance to the neighbours

        lap[i, i] = -(4 / hi) * invhi  # Laplacian of vertex itself
        lap[i, k] = (4 / (hi * ni)) * 1. / edge[i, k]  # Laplacian of direct neighbours

        # Laplacian is zero for all indirect neighbours
        # See Oostendorp, Oosterom & Huiskamp (1989, pp. 334-335)

    edge = csr_array(edge)
    lap = csr_array(lap)

    print(f'done.')
    
    return lap, edge

def mesh_laplacian_interp(lap, index):
    '''
     MESH_LAPLACIAN_INTERP: Computes the zero Laplacian interpolation matrix
 
    Usage:   mesh_laplacian_interp(lap, index) -> [int, keepindex, repindex]
    
    This function calculates an interpolation matrix that provides
    the coefficients for the calculation of potential values at all
    unknown vertices of a mesh, given known potential values at
    a subset of the mesh vertices (at 'index').  The interpolation
    solution is constrained by a minimal norm of the Laplacian
    of the mesh.  See the reference below for details.
    
    'lap' is the laplacian matrix for the full mesh (see mesh_laplacian)
    'int' is the matrix which interpolates from the points in 'index'
    to the full mesh.  'index' is a row vector of indices into a 
    subset of the vertices used to calculate 'lap'.  This subset 
    is where the electric potential is known and usually corresponds 
    to the given electrode vertices, eg:
    
    index = dsearchn(scalpvert,elecvert)';
    
    If 'index' contains repeated indices, only the unique indices 
    are useful.  The 'keepindex' array can be used to select these.
    The 'repindex' array is the repeated indices.
    
    Interpolations can be done using matrix 'int', eg:
    
    [int, keepindex, repindex] = mesh_laplacian_interp(lap,index);
    if isempty(repindex),
    Vint = int * Vknown;
    else
    Vint = int * Vknown(keepindex);
    end
    
    This implements interpolation method B (p. 336) of 
    Oostendorp T, Oosterom A & Huiskamp G (1989),
    Interpolation on a triangulated 3D surface.
    Journal of Computational Physics, 80: 331-343.


    Licence:  GNU GPL, no implied or express warranties
    History:  (c) 04/2002 Robert Oostenveld
            - agreed to release 'lapint' under GNU GPL
            04/2002, Darren.Weber@flinders.edu.au
            - introduced check for index replications and
                adjusted calculations/output accordingly
            - converted lap to sparse matrix and solution
                of interpolation matrix with \ operator
            - accepts sparse lap input and returns sparse int

    '''


    if lap.shape[0] != lap.shape[1]:
        raise ValueError("MESH_LAPLACIAN_INTERP: lap matrix is not square")

    if issparse(lap):
        lap = lap.toarray()


    known_index = np.sort(np.unique(index))
    k = len(known_index)
    if len(index) != k: print("THE INPUT INDEX HAS REPEATED VALUES!")
           
    n = lap.shape[0]

    print(f'MESH_LAPLACIAN_INTERP: Calc Interpolation matrix for {k} to {n} vertices... ', end="")
    
    # find 'unknown' indices of lap matrix
    unknown_index = np.setdiff1d(np.arange(n), known_index)

    # reshuffle rows & columns of lap matrix
    lapi = np.concatenate((known_index, unknown_index))
    lap = lap[lapi, :]  # rows
    lap = lap[:, lapi]  # columns

    # Segregate known/unknown portions of lap
    L11 = lap[:k, :k]
    L12 = lap[:k, k:n]
    L21 = lap[k:n, :k]
    L22 = lap[k:n, k:n]

    # Convert to sparse for quicker computation
    A = csr_array(np.concatenate((L12, L22)))
    B = csr_array(np.concatenate((L11, L21)))

    # Perform sparse matrix division
    int_matrix = np.matmul(np.linalg.pinv(-A.toarray()), B.toarray())

    # Convert result back to full matrix
    # int_matrix = int_matrix.toarray()

    # append the interpolating piece to the identity matrix
    # these take care of the known potentials
    int_matrix = np.block([[np.eye(k)], [int_matrix]])

    # reshuffle the columns of the interpolating matrix
    order = np.argsort(known_index)
    int_matrix = int_matrix[:, order]

    # reshuffle the rows of the interpolating matrix
    order = np.argsort(np.concatenate((known_index, unknown_index)))
    int_matrix = int_matrix[order, :]

    print('done.')
    
    return int_matrix, known_index

def mesh_blur(signals, faces, it=1):
    connectivity = get_connectivity(faces) 
    for _ in range(it):
        outsignals = np.zeros_like(signals)
        
        for i in range(signals.shape[0]):
            outsignals[i] = signals[connectivity[i]].mean(axis=0)
            
        signals = outsignals
    return signals