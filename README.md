# spatialHMRF

Spatial domains with the hidden Markov random field (HMRF) model of
[Zhu et al. 2018](https://doi.org/10.1038/nbt.4260): Gaussian emissions on
z-scored spatial genes, a Potts prior over the spatial network, k-means start,
and EM with iterated conditional modes over the colour classes of the network.

It reimplements the model of the `smfishHmrf` package (used by Giotto's
`doHMRF_V2()`) in plain R with sparse matrix products, and gives the same
results on the same inputs. It depends only on Matrix and igraph.

```r
fit <- spatialHMRF::hmrf(y, network, k = 8, beta = 2)
fit$class   # domain per cell
fit$prob    # posterior probabilities
```

In Giotto, use `runHMRF()` or `findSpatialDomains(method = "hmrf")`.

**If you use it, please cite** Zhu Q, Shah S, Dries R, Cai L, Yuan GC.
Identification of spatially associated subpopulations by combining scRNAseq
and sequential fluorescence in situ hybridization data. *Nature Biotechnology*
36, 1183–1190 (2018).
