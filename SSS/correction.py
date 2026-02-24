from os.path import join
from glob import glob

import numpy as np
import pandas as pd

from scipy.stats import ttest_1samp

from tqdm import tqdm

from SSS import util as su
from SSS import deal_spm
from SSS import image as simage
from SSS import stat as sstat

import nibabel as nb

def get_geometry_from_fs32k(fname):
    """
    Load a GIFTI (.gii) surface file and extract its geometric information.

    Parameters
    ----------
    fname : str
        Path to the GIFTI file to be loaded.

    Returns
    -------
    coords : ndarray of shape (N, 3)
        3D coordinates of each vertex in the surface mesh.
    faces : ndarray of shape (M, 3)
        Triangle indices defining the mesh connectivity. Each row contains
        the vertex indices that form one triangular face.
    """
    gii = nb.load(fname)
    coords = gii.darrays[0].data
    faces = gii.darrays[1].data.astype(int)

    return coords, faces

def get_adjacenct_matrix_fs32k(fname):
    """
    Construct an adjacency list for a fsaverage32k surface mesh.

    Parameters
    ----------
    fname : str
        Path to the GIFTI (.gii) surface file containing geometry information.

    Returns
    -------
    adj : list of lists
        Adjacency list where adj[i] contains all unique vertices that share
        a face with vertex i. Each inner list corresponds to one vertex's
        neighbors in the mesh.
    """
    coords, faces = get_geometry_from_fs32k(fname)
    V = coords.shape[0]
    tmp = [[] for _ in range(V)]
    for f in faces:
        a, b, c = map(int, f)
        tmp[a].extend([b,c])
        tmp[b].extend([a,c])
        tmp[c].extend([a,b])
    adj = [list(set(inner)) for inner in tmp]

    return adj

def find_clusters(supra, adj):
    """
    Identify connected clusters of vertices where `supra` is True,
    using depth-first search (DFS) on a surface adjacency structure.

    Parameters
    ----------
    supra : array-like of bool, shape (V,)
        Boolean mask indicating which vertices are "active" or suprathreshold.
        Only vertices with supra[v] == True are considered for clustering.

    adj : list of lists
        Adjacency list where adj[v] contains the neighboring vertices of v.

    Returns
    -------
    clusters : list of lists
        A list of clusters, where each cluster is a list of vertex indices
        forming one connected component among the True entries of `supra`.
    """
    V = len(supra)
    visited = np.zeros(V, dtype=bool)
    clusters = []

    for v in np.where(supra)[0]:
        v = int(v)
        if visited[v]:
            continue
        stack = [v]
        cluster = []
        while stack:
            u = stack.pop()
            if visited[u] or not supra[u]:
                continue
            visited[u] = True
            cluster.append(u)
            for w in adj[u]:
                if not visited[w] and supra[w]:
                    stack.append(w)
        if len(cluster) > 0:
            clusters.append(cluster)

    return clusters

    # # example: 0-1-2 | 3-4-5 | 6-7 | 8
    # adj = { 0:[1], 1:[0, 2], 2:[1], 3:[4], 4:[3,5], 5:[4], 6:[7], 7:[6], 8:[] }
    # # vertex 2 and 7 < thresh
    # supra = np.array([True, True, False, True, True, True, True, False, True])
    # # 함수 실행
    # clusters = find_clusters(supra, adj)
    # print("Detected clusters:", clusters)

def compute_null_max_cluster_mass(X, thresh, adj, n_perm=5000, alternative='greater', axis=0, popmean=0):
    """
    Generate a null distribution of maximum cluster masses using sign-flip
    permutation testing on surface-based data.

    Parameters
    ----------
    X : ndarray, shape (N_subjects, V)
        Input data matrix where each row is a subject and each column is a vertex.

    thresh : float
        Threshold applied to the t-statistic to define suprathreshold vertices.

    adj : list of lists
        Surface adjacency list. adj[v] contains neighbors of vertex v.

    n_perm : int, default=5000
        Number of permutation iterations.

    alternative : {'greater', 'less', 'two-sided'}, default='greater'
        Alternative hypothesis for the one-sample t-test.

    axis : int, default=0
        Axis along which the t-test is computed (subjects axis).

    popmean : float, default=0
        Population mean for the one-sample t-test.

    Returns
    -------
    clusters : list of lists
        Clusters detected in the *final* permutation iteration. Mainly useful
        for debugging; not typically used for inference.
    null_dist : ndarray of shape (n_perm,)
        Null distribution of maximum cluster masses across permutations.
    """
    max_cluster_masses = []
    for _ in tqdm(range(n_perm)):
        signs = np.random.choice([1, -1], size=(12, 1))
        X_perm = X * signs

        t_perm, _ = ttest_1samp(X_perm, popmean=popmean, axis=axis, alternative=alternative)
    
        supra_perm = t_perm > thresh

        clusters = find_clusters(supra_perm, adj)

        if len(clusters) > 0: 
            max_mass = max([np.max(t_perm[c]) for c in clusters]) 
        else: 
            max_mass = 0 
            
        max_cluster_masses.append(max_mass)
        
    return clusters, np.array(max_cluster_masses)
    
def compute_cluster_pvals(tmap, thresh, null_dist, adj):
    """
    Compute p-values for suprathreshold clusters in a t-statistic map
    using a permutation-derived null distribution of cluster masses.

    Parameters
    ----------
    tmap : ndarray of shape (V,)
        Vertex-wise t-statistic map.

    thresh : float
        Threshold applied to the t-statistic to define suprathreshold vertices.

    null_dist : ndarray of shape (n_perm,)
        Null distribution of maximum cluster masses obtained from permutation testing.

    adj : list of lists
        Surface adjacency list. adj[v] contains neighbors of vertex v.

    Returns
    -------
    clusters : list of lists
        List of detected clusters, where each cluster is a list of vertex indices.

    pvals : ndarray of shape (n_clusters,)
        Cluster-wise p-values computed as the proportion of null masses
        greater than or equal to each real cluster's mass.
    """
    supra = tmap > thresh
    clusters = find_clusters(supra, adj)
    real_masses = np.array([np.sum(tmap[c]) for c in clusters])
    real_masses = np.array([sum(tmap[c]) for c in clusters])
    pvals = np.array([np.mean(null_dist >= mass) for mass in real_masses])

    return clusters, pvals

def make_cluster_corrected_mask(tmap, clusters, pvals, alpha=0.5, fvalue=np.nan):
    """
    Create a cluster-corrected mask by retaining only clusters whose p-values
    survive a specified significance threshold.

    Parameters
    ----------
    tmap : ndarray of shape (V,)
        Vertex-wise t-statistic map.

    clusters : list of lists
        List of clusters, where each cluster is a list of vertex indices.

    pvals : ndarray of shape (n_clusters,)
        Cluster-wise p-values corresponding to the clusters.

    alpha : float, default=0.5
        Significance threshold. Clusters with p < alpha are retained.

    fvalue : float, default=np.nan
        Fill value assigned to non-significant vertices.

    Returns
    -------
    mask : ndarray of shape (V,)
        Cluster-corrected mask. Significant clusters contain 1, and all
        non-significant vertices are filled with `fvalue`.
    """
    masked = np.ones_like(tmap) * np.nan
    for cluster, p in zip(clusters, pvals):
        if p < alpha:
            masked[cluster] = tmap[cluster]

    mask = np.where(np.isnan(masked), fvalue, 1)

    return mask