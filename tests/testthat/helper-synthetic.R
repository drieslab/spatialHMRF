# Two-domain tissue with a symmetric kNN network.
sim_tissue <- function(n = 300, seed = 1, k_nn = 5) {
    set.seed(seed)
    xy <- matrix(runif(2 * n), ncol = 2) * 10
    dom <- 1L + (xy[, 1] > 5)
    y <- matrix(rnorm(n * 6), n) + c(-1.5, 1.5)[dom]
    ids <- paste0("c", seq_len(n))
    rownames(y) <- ids
    d <- as.matrix(stats::dist(xy))
    diag(d) <- Inf
    knn <- t(apply(d, 1, order))[, seq_len(k_nn)]
    e <- unique(t(apply(cbind(rep(seq_len(n), k_nn), as.vector(knn)), 1, sort)))
    list(y = y, dom = dom, net = data.frame(from = ids[e[, 1]], to = ids[e[, 2]]))
}
ari <- function(a, b) {
    tab <- table(a, b)
    comb <- function(x) sum(choose(x, 2))
    n <- sum(tab)
    exp <- comb(rowSums(tab)) * comb(colSums(tab)) / choose(n, 2)
    (comb(tab) - exp) / ((comb(rowSums(tab)) + comb(colSums(tab))) / 2 - exp)
}
