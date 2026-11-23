# Equivalence and wall time of spatialHMRF vs smfishHmrf on simulated tissue (3 domains,
# 20 genes, symmetric 6-NN network), same inputs, beta = 2, k = 3.
suppressPackageStartupMessages({ library(spatialHMRF); library(smfishHmrf) })
sim <- function(n, seed = 1) {
    set.seed(seed)
    xy <- matrix(runif(2 * n), ncol = 2) * sqrt(n)
    dom <- 1L + (xy[, 1] > sqrt(n) / 3) + (xy[, 1] > 2 * sqrt(n) / 3)
    mu <- matrix(rnorm(3 * 20, sd = 0.6), 3)
    y <- mu[dom, ] + matrix(rnorm(n * 20), n)
    ids <- paste0("c", seq_len(n)); rownames(y) <- ids
    knn <- dbscan::kNN(xy, k = 6)$id
    e <- unique(t(apply(cbind(rep(seq_len(n), 6), as.vector(knn)), 1, sort)))
    list(y = y, dom = dom, net = data.frame(from = ids[e[, 1]], to = ids[e[, 2]]))
}
res <- list()
for (n in c(2000, 10000, 50000, 100000)) {
    s <- sim(n); ids <- rownames(s$y)
    adj <- hmrf_adjacency(s$net, ids); blocks <- hmrf_blocks(adj)
    init <- hmrf_init(s$y, 3, nstart = 10)
    t0 <- proc.time()[[3]]; nat <- hmrf(s$y, s$net, k = 3, beta = 2, zscore = "none", init = init, blocks = blocks); t_nat <- proc.time()[[3]] - t0
    t_ref <- NA; same <- NA
    if (n <= 10000) {
        nb <- split(c(s$net$to, s$net$from), c(s$net$from, s$net$to))[ids]
        nl <- lapply(nb, function(v) sort(match(v, ids))); numnei <- lengths(nl)
        nei <- matrix(-1L, n, max(numnei)); for (i in seq_len(n)) nei[i, seq_len(numnei[i])] <- nl[[i]]
        t0 <- proc.time()[[3]]
        ref <- smfishHmrf.hmrfem.multi(y = s$y, neighbors = nei, beta = 2, numnei = as.integer(numnei),
            blocks = blocks, mu = init$mu, sigma = init$sigma, verbose = FALSE, err = 1e-7,
            maxit = 50, dampFactor = init$damp)
        t_ref <- proc.time()[[3]] - t0
        same <- identical(as.integer(ref$class), as.integer(nat$class))
    }
    ari <- aricode::ARI(nat$class, s$dom)
    cat(sprintf("n=%6d  native %.1fs  smfishHmrf %s  identical=%s  ARI(truth)=%.3f  iters=%d\n",
        n, t_nat, if (is.na(t_ref)) "skipped" else sprintf("%.1fs", t_ref), same, ari, nat$iterations))
    res[[length(res) + 1]] <- data.frame(n = n, t_native = t_nat, t_smfishHmrf = t_ref, identical = same, ari = ari)
}
dir.create("results", showWarnings = FALSE)
write.csv(do.call(rbind, res), "results/scaling.csv", row.names = FALSE)
