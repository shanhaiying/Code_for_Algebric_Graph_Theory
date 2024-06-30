
# -*- coding: utf-8 -*-

#########################################################################################
#      This script enhances the NetworkX library with advanced spectral graph theory    #
#      functions, enabling computations related to various graph matrices and their     #
#      eigenvalues.                                                                     #
#                                                                                       #
#       Copyright (C) 2014~2024 HaiyingShan <shan_haiying@tongji.edu.cn>                #
#                             Last Modified: 2024-6-30                                #
#########################################################################################

import networkx as nx
import numpy as np
import sympy as sp
from numpy.polynomial import Polynomial

def is_integer_matrix(M):
    return np.all(np.mod(M, 1) == 0)

def degree_matrix(G):
    degrees = np.array([d for n, d in G.degree()])
    D = np.diag(degrees)
    return D
nx.Graph.degree_matrix = degree_matrix

def signless_laplacian_matrix(G):
    D = degree_matrix(G)
    A = nx.adjacency_matrix(G).todense()
    return D + A

def alpha_matrix(G, alpha):
    A = nx.adjacency_matrix(G).todense()
    D = degree_matrix(G)
    return alpha * D + (1 - alpha) * A

def distance_alpha_matrix(G, alpha):
    D = nx.floyd_warshall_numpy(G)
    T = np.diag(np.sum(D, axis=0))
    return alpha * T + (1 - alpha) * D

def abc_matrix(G):
    n = G.order()
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if G.has_edge(i, j):
                di = G.degree(i)
                dj = G.degree(j)
                M[i, j] = np.sqrt((di + dj - 2) / (di * dj))
    return M

def dist_lap(G):
    n = G.order()
    ones = np.ones((n, 1))
    dg = nx.floyd_warshall_numpy(G)
    trg = dg @ ones
    DL = np.diag(np.asarray(trg).flatten()) - dg
    return DL

def hermitian_adjacency_matrix(G):
    n = G.order()
    I = 1j
    A = nx.adjacency_matrix(G).todense()
    H = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            if G.has_edge(i, j):
                H[i, j] = I if not G.has_edge(j, i) else 1
                if not G.has_edge(j, i):
                    H[j, i] = -I
    return H

# Extending NetworkX Graph and DiGraph classes
nx.Graph.signless_laplacian_matrix = signless_laplacian_matrix
nx.Graph.alpha_matrix = alpha_matrix
nx.DiGraph.alpha_matrix = alpha_matrix
nx.Graph.distance_alpha_matrix = distance_alpha_matrix
nx.Graph.abc_matrix = abc_matrix
nx.Graph.dist_lap = dist_lap
nx.DiGraph.hermitian_adjacency_matrix = hermitian_adjacency_matrix

def eigenpair(M, spe=False, order=0):
    """
    Calculate the k-th largest eigenvalue and corresponding eigenvector of a real symmetric matrix.
    
    Parameters:
    M (numpy.ndarray): Real symmetric matrix
    spe (bool): Whether to return all eigenvalues
    order (int or str): Order of the eigenvalue to return, 0 for largest, 1 for second largest, etc.;
                        "all" to return all eigenvalues and eigenvectors
    
    Returns:
    tuple or list: k-th largest eigenvalue and corresponding eigenvector, or all eigenvalues and eigenvectors
    """
    if not np.allclose(M, M.T):
        return "A symmetric matrix is required."
    
    w, v = np.linalg.eigh(M)
    v = np.around(v, decimals=12)
    w = np.around(w, decimals=12)
    
    idx = np.argsort(w)[::-1]
    w = w[idx]
    v = v[:, idx]
    
    if spe:
        return w
    if order == "all":
        return w, v
    else:
        rec = v[:, order]
        if rec[0] != 0:
            rec = np.sign(rec[0]) * rec
        return w[order], rec

def eigenpairs(self, type="L", order=0, spe=False):
    if type == "A":
        M = nx.adjacency_matrix(self).todense()
    elif type == "L":
        M = nx.laplacian_matrix(self).todense()
    elif type == "N":
        M = nx.normalized_laplacian_matrix(self).todense()
    elif type == "Q":
        M = self.signless_laplacian_matrix()
    else:
        raise ValueError(f"Unknown matrix type: {type}")
    return eigenpair(M, spe=spe, order=order)

nx.Graph.eigenpairs = eigenpairs

def charpolys_for_graph(G, type="L", alpha=0.5):
    if isinstance(type, (int, float)):
        M = alpha_matrix(G, type)
    elif type == "A":
        M = nx.adjacency_matrix(G).todense()
    elif type == "L":
        M = nx.laplacian_matrix(G).todense()
    elif type == "N":
        M = nx.normalized_laplacian_matrix(G).todense()
    elif type == "Q":
        M = signless_laplacian_matrix(G)
    elif type == "DA":
        M = distance_alpha_matrix(G, alpha)
    elif type == "ABC":
        M = abc_matrix(G)
    elif type == "DL":
        M = dist_lap(G)
    elif type == "H":
        M = hermitian_adjacency_matrix(G)
    else:
        raise ValueError(f"Unknown matrix type: {type}")
    charpoly_coeffs = np.poly(M)
    if is_integer_matrix(M):
        charpoly_coeffs = charpoly_coeffs.astype(int)
    charpoly = Polynomial(charpoly_coeffs[::-1])
    return charpoly
    
nx.Graph.charpolys_for_graph = charpolys_for_graph
nx.Graph.charpolys = charpolys_for_graph
nx.DiGraph.charpolys = charpolys_for_graph
