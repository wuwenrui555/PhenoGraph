from typing import Union, Optional, Type

import leidenalg
import igraph as ig
import numpy as np
from leidenalg.VertexPartition import MutableVertexPartition
from scipy import sparse as sp
from scipy.sparse.base import spmatrix

from phenograph.core import (
    gaussian_kernel,
    parallel_jaccard_kernel,
    jaccard_kernel,
    find_neighbors,
    neighbor_graph,
    graph2binary,
    runlouvain,
)
import time
import re
import os
import uuid


def sort_by_size(clusters, min_size):
    """
    Relabel clustering in order of descending cluster size.
    New labels are consecutive integers beginning at 0
    Clusters that are smaller than min_size are assigned to -1
    :param clusters:
    :param min_size:
    :return: relabeled
    """
    relabeled = np.zeros(clusters.shape, dtype=np.int)
    sizes = [sum(clusters == x) for x in np.unique(clusters)]
    o = np.argsort(sizes)[::-1]
    for i, c in enumerate(o):
        if sizes[c] > min_size:
            relabeled[clusters == c] = i
        else:
            relabeled[clusters == c] = -1
    return relabeled


def cluster(
    data: Union[np.ndarray, spmatrix],
    clustering_algo: Union["louvain", "leiden"] = "louvain",
    k: int = 30,
    directed: bool = False,
    prune: bool = False,
    min_cluster_size: int = 10,
    jaccard: bool = True,
    primary_metric: Union[
        "euclidean", "manhattan", "correlation", "cosine"
    ] = "euclidean",
    n_jobs: int = -1,
    q_tol: float = 1e-3,
    louvain_time_limit: int = 2000,
    nn_method: Union["kdtree", "brute"] = "kdtree",
    partition_type: Optional[Type[MutableVertexPartition]] = None,
    resolution_parameter: float = 1,
    use_weights: bool = True,
    seed: int = 0,
):
    """
    PhenoGraph clustering

    :param data: Numpy ndarray of data to cluster, or sparse matrix of k-nearest
        neighbor graph.
        If ndarray, n-by-d array of n cells in d dimensions
        If sparse matrix, n-by-n adjacency matrix
    :param clustering_algo: Optional `'louvain'`, or `'leiden'`. Any other value will
        return only graph object.
    :param k: Number of nearest neighbors to use in first step of graph construction
    :param directed: Whether to use a symmetric (default) or asymmetric ("directed")
        graph
        The graph construction process produces a directed graph, which is symmetrized
        by one of two methods (see below)
    :param prune: Whether to symmetrize by taking the average (prune=False) or product
        (prune=True) between the graph
        and its transpose
    :param min_cluster_size: Cells that end up in a cluster smaller than
        min_cluster_size are considered outliers and are assigned to -1 in the cluster
        labels
    :param jaccard: If True, use Jaccard metric between k-neighborhoods to build graph.
        If False, use a Gaussian kernel.
    :param primary_metric: Distance metric to define nearest neighbors.
        Options include: {'euclidean', 'manhattan', 'correlation', 'cosine'}
        Note that performance will be slower for correlation and cosine.
    :param n_jobs: Nearest Neighbors and Jaccard coefficients will be computed in
        parallel using n_jobs. If n_jobs=-1, the number of jobs is determined
        automatically.  This parameter is passed to Leiden as `n_iterations`, the number
        of iterations to run the Leiden algorithm. If the number of iterations is
        negative, the Leiden algorithm is run until an iteration in which there was no
        improvement.
    :param q_tol: Tolerance (i.e., precision) for monitoring modularity optimization
    :param louvain_time_limit: Maximum number of seconds to run modularity optimization.
        If exceeded the best result so far is returned
    :param nn_method: Whether to use brute force or kdtree for nearest neighbor search.
        For very large high-dimensional data sets, brute force (with parallel
        computation) performs faster than kdtree.
    :param partition_type: Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    :param resolution_parameter: A parameter value controlling the coarseness of the
        clustering in Leiden. Higher values lead to more clusters. Set to `None` if
        overriding `partition_type` to one that doesn’t accept a `resolution_parameter`.
    :param use_weights: Use vertices in the Leiden computation.
    :param seed: Leiden initialization of the optimization

    :return communities: numpy integer array of community assignments for each row in data
    :return graph: numpy sparse array of the graph that was used for clustering
    :return Q: the modularity score for communities on graph
    """

    # NB if prune=True, graph must be undirected, and the prune setting takes precedence
    if prune:
        print("Setting directed=False because prune=True")
        directed = False

    if n_jobs == 1:
        kernel = jaccard_kernel
    else:
        kernel = parallel_jaccard_kernel
    kernelargs = {}

    # Start timer
    tic = time.time()
    # Go!
    if isinstance(data, sp.spmatrix) and data.shape[0] == data.shape[1]:
        print(
            "Using neighbor information from provided graph, "
            "rather than computing neighbors directly",
            flush=True,
        )
        lilmatrix = data.tolil()
        d = np.vstack(lilmatrix.data).astype("float32")  # distances
        idx = np.vstack(lilmatrix.rows).astype("int32")  # neighbor indices by row
        del lilmatrix
        assert idx.shape[0] == data.shape[0]
        k = idx.shape[1]
    else:
        d, idx = find_neighbors(
            data, k=k, metric=primary_metric, method=nn_method, n_jobs=n_jobs
        )
        print("Neighbors computed in {} seconds".format(time.time() - tic), flush=True)

    subtic = time.time()
    kernelargs["idx"] = idx
    # if not using jaccard kernel, use gaussian
    if not jaccard:
        kernelargs["d"] = d
        kernelargs["sigma"] = 1.0
        kernel = gaussian_kernel
        graph = neighbor_graph(kernel, kernelargs)
        print(
            "Gaussian kernel graph constructed in {} seconds".format(
                time.time() - subtic
            ),
            flush=True,
        )
    else:
        del d
        graph = neighbor_graph(kernel, kernelargs)
        print(
            "Jaccard graph constructed in {} seconds".format(time.time() - subtic),
            flush=True,
        )
    if not directed:
        if not prune:
            # symmetrize graph by averaging with transpose
            sg = (graph + graph.transpose()).multiply(0.5)
        else:
            # symmetrize graph by multiplying with transpose
            sg = graph.multiply(graph.transpose())
        # retain lower triangle (for efficiency)
        graph = sp.tril(sg, -1)

    # choose between Louvain or Leiden algorithm
    communities, Q = "", ""
    if clustering_algo == "louvain":
        # write to file with unique id
        uid = uuid.uuid1().hex
        graph2binary(uid, graph)
        communities, Q = runlouvain(uid, tol=q_tol, time_limit=louvain_time_limit)
        print("PhenoGraph complete in {} seconds".format(time.time() - tic), flush=True)
        communities = sort_by_size(communities, min_cluster_size)
        # clean up
        for f in os.listdir():
            if re.search(uid, f):
                os.remove(f)

    elif clustering_algo == "leiden":
        # convert resulting graph from scipy.sparse.coo.coo_matrix to Graph object
        # get indices of vertices
        edgelist = np.vstack(graph.nonzero()).T.tolist()
        g = ig.Graph(max(graph.shape), edgelist, directed=directed)
        # set vertices as weights
        g.es["weights"] = graph.toarray()[graph.nonzero()]

        kargs = dict()
        if not partition_type:
            partition_type = leidenalg.RBConfigurationVertexPartition
        if resolution_parameter:
            kargs["resolution_parameter"] = resolution_parameter
        if use_weights:
            kargs["weights"] = np.array(g.es["weights"]).astype("float64")
        kargs["n_iterations"] = n_jobs
        kargs["seed"] = seed

        print("Running Leiden optimization", flush=True)
        tic_ = time.time()
        communities = leidenalg.find_partition(
            g, partition_type=partition_type, **kargs,
        )
        Q = communities.q
        print(
            "Leiden completed in {} seconds".format(time.time() - tic_), flush=True,
        )
        communities = np.asarray(communities.membership)
        print("PhenoGraph complete in {} seconds".format(time.time() - tic), flush=True)
        communities = sort_by_size(communities, min_cluster_size)

    else:
        # return only graph object
        pass

    return communities, graph, Q
