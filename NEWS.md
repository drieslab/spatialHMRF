# spatialHMRF 0.1.1

* `zscore = "feature"` scales each feature across cells, for a few features
  on different scales such as module scores.
* Fits survive degenerate domains: one-cell start clusters, emptied domains,
  and covariances that lose positive definiteness (the damping grows).

# spatialHMRF 0.1.0

* First version: `hmrf()`, `hmrf_init()`, `hmrf_adjacency()`, `hmrf_blocks()`.
* k-means start: `nstart` restarts on a `subsample` of cells, the best 10
  refined on all cells.
* EM stops when fewer than `label_tol` of the cells change domain;
  densities through Cholesky factors.
