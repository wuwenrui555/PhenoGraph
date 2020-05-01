PhenoGraph for Python3
======================

[PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) is a clustering method designed for 
high-dimensional single-cell data. It works by creating a graph ("network") representing phenotypic similarities 
between cells and then identifying communities in this graph. 

This implementation is written in Python3 and depends only on `scikit-learn (>= 0.17)` and its dependencies.  

This software package includes compiled binaries that run community detection based on C++ code written by 
E. Lefebvre and J.-L. Guillaume in 2008 (["Louvain method"](https://sites.google.com/site/findcommunities/)). The code
has been altered to interface more efficiently with the Python code here. It should work on reasonably current Linux, 
Mac and Windows machines.

To install PhenoGraph, simply run the setup script:

    pip install PhenoGraph


Expected use is within a script or interactive kernel running Python `3.x`. Data are expected to be passed as a `numpy.ndarray`.
When applicable, the code uses CPU multicore parallelism via `multiprocessing`. 
 
To run basic clustering:

    import phenograph
    communities, graph, Q = phenograph.cluster(data)

For a dataset of *N* rows, `communities` will be a length *N* vector of integers specifying a community assignment for each row
in the data. Any rows assigned `-1` were identified as *outliers* and should not be considered as a member of any community.
`graph` is a *N* x *N* `scipy.sparse` matrix representing the weighted graph used for community detection. 
`Q` is the modularity score for `communities` as applied to `graph`.

If you use PhenoGraph in work you publish, please cite our publication:

    @article{Levine_PhenoGraph_2015,
      doi = {10.1016/j.cell.2015.05.047},
      url = {http://dx.doi.org/10.1016/j.cell.2015.05.047},
      year  = {2015},
      month = {jul},
      publisher = {Elsevier {BV}},
      volume = {162},
      number = {1},
      pages = {184--197},
      author = {Jacob H. Levine and Erin F. Simonds and Sean C. Bendall and Kara L. Davis and El-ad D. Amir and Michelle D. Tadmor and Oren Litvin and Harris G. Fienberg and Astraea Jager and Eli R. Zunder and Rachel Finck and Amanda L. Gedman and Ina Radtke and James R. Downing and Dana Pe'er and Garry P. Nolan},
      title = {Data-Driven Phenotypic Dissection of {AML} Reveals Progenitor-like Cells that Correlate with Prognosis},
      journal = {Cell}
    }

Release Notes
-------------

### Version 1.5.5
 * Exposed parameter `n_iterations` for Leiden, along with minor fixes to manage sorting parallelism, and updated documentation of the clustering and sorting methods.

### Version 1.5.4

 * Faster and more efficient sorting by size of clusters, for large nearest neighbours graph, implementing multiprocessing and faster methods for sorting.
 
### Version 1.5.3

 * Phenograph supports now [**Leiden**](https://www.nature.com/articles/s41598-019-41695-z) algorithm for community detection.
 The new feature can be called from `phenograph.cluster`, by choosing `leiden` as the clustering algorithm. 
 
### Version 1.5.2

 * Include simple parallel implementation of brute force nearest neighbors search using scipy's `cdist` and `multiprocessing`. This may be more efficient than `kdtree` on very large high-dimensional data sets
 and avoids memory issues that arise in `sklearn`'s implementation.
 * Refactor `parallel_jaccard_kernel` to remove unnecessary use of `ctypes` and `multiprocessing.Array`.

### Version 1.5.1

 * Make `louvain_time_limit` a parameter to `phenograph.cluster`.

### Version 1.5

 * `phenograph.cluster` can now take as input a square sparse matrix, which will be interpreted as a k-nearest neighbor graph. 
 Note that this graph _must_ have uniform degree (i.e. the same value of k at every point).
 * The default `time_limit` for Louvain iterations has been increased to a more generous 2000 seconds (~half hour).
  
### Version 1.4.1

 * After observing inconsistent behavior of sklearn.NearestNeighbors with respect to inclusion of self-neighbors,
 the code now checks that self-neighbors have been included before deleting those entries.
 
### Version 1.4

 * The dependence on IPython and/or ipyparallel has been removed. Instead the native `multiprocessing` package is used.
 * Multiple CPUs are used by default for computation of nearest neighbors and Jaccard graph.
 
### Version 1.3

 * Proper support for Linux.

---
Troubleshooting
---------------

### Notebook freezes after several attempts of running PhenoGraph using Jypyter Notebook
 
* Running `PhenoGraph` from a Jupyter Notebook repeatedly on macOS Catalina, but not Mojave,  using Python 3.7.6, causes a hang and the notebook becomes unresponsive, even for a basic matrix of nearest neighbors. However, this issue was not reproducible in command line using `Python` interpreter in both Catalina and Mojave platforms, without using Jupyter Notebook.
  
  It was found that attempting to plot principal components using 
    
    ```
    :func:`~matplotlib.pyplot.scatter`
    ```
    
  in Jupyter Notebook causes a freeze, and `PhenoGraph` becomes unresponsive unless the kernel is restarted. When removing this line of code, everything goes back to normal and the Jupyter Notebook stopes crashing with repeated runs of `PhenoGraph`. 
  
### Architecture related error

* When attempting to process very large nearest neighbours graph, _e.g._ a 2000000 `x` 2000000 kNN graph matrix with 300 nearest neighbours, a `struct.error()` is raised: 
  
    ```python
    struct.error: 'i' format requires -2147483648 <= number <= 2147483647
    ```
  
  This issue was reported on [stackoverflow](https://stackoverflow.com/questions/47776486/python-struct-error-i-format-requires-2147483648-number-2147483647) and it's related to the multiprocessing while building the Jaccard object.
  
  The `struct.error()` has been fixed in python >= 3.8.0.
    
### `leidenalg` inside conda environment

* When using `PhenoGraph` inside a conda environment `leiden` takes longer to complete for larger samples compared to the system Python. 

