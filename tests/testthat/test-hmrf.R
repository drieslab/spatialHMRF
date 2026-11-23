test_that("hmrf finds the simulated domains", {
    s <- sim_tissue()
    fit <- hmrf(s$y, s$net, k = 2, beta = 2, zscore = "none", nstart = 10)
    expect_length(fit$class, nrow(s$y))
    expect_equal(dim(fit$prob), c(nrow(s$y), 2L))
    expect_equal(unname(rowSums(fit$prob)), rep(1, nrow(s$y)))
    expect_equal(rownames(fit$prob), rownames(s$y))
    expect_gt(ari(fit$class, s$dom), 0.9)
})

test_that("blocks are independent sets covering every cell once", {
    s <- sim_tissue()
    adj <- hmrf_adjacency(s$net, rownames(s$y))
    b <- hmrf_blocks(adj)
    expect_setequal(unlist(b), seq_len(nrow(s$y)))
    for (bl in b) expect_equal(sum(adj[bl, bl]), 0)
})

test_that("a start and blocks can be reused across beta", {
    s <- sim_tissue()
    f1 <- hmrf(s$y, s$net, k = 2, beta = 0, zscore = "none", nstart = 5)
    f2 <- hmrf(s$y, s$net, k = 2, beta = 4, zscore = "none",
        init = f1$init, blocks = f1$blocks)
    expect_identical(f2$init, f1$init)
    expect_gte(f2$likelihood, -Inf)
})

test_that("cells without neighbours are handled", {
    s <- sim_tissue()
    net <- s$net[s$net$from != "c1" & s$net$to != "c1", ]
    fit <- hmrf(s$y, net, k = 2, beta = 2, zscore = "none", nstart = 5)
    expect_false(anyNA(fit$class))
})

test_that("matches smfishHmrf on identical inputs", {
    skip_if_not_installed("smfishHmrf")
    s <- sim_tissue()
    ids <- rownames(s$y)
    adj <- hmrf_adjacency(s$net, ids)
    blocks <- hmrf_blocks(adj)
    init <- hmrf_init(s$y, 2, nstart = 10)
    nb <- split(c(s$net$to, s$net$from), c(s$net$from, s$net$to))[ids]
    nl <- lapply(nb, function(v) sort(match(v, ids)))
    numnei <- lengths(nl)
    nei <- matrix(-1L, length(ids), max(numnei))
    for (i in seq_along(ids)) nei[i, seq_len(numnei[i])] <- nl[[i]]
    ref <- smfishHmrf::smfishHmrf.hmrfem.multi(y = s$y, neighbors = nei,
        beta = 2, numnei = as.integer(numnei), blocks = blocks,
        mu = init$mu, sigma = init$sigma, verbose = FALSE, err = 1e-7,
        maxit = 50, dampFactor = init$damp)
    fit <- hmrf(s$y, s$net, k = 2, beta = 2, zscore = "none",
        init = init, blocks = blocks, label_tol = 0)
    expect_identical(as.integer(fit$class), as.integer(ref$class))
    expect_equal(fit$prob, ref$prob, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(fit$likelihood, ref$likelihood, tolerance = 1e-10)
})

test_that("label_tol stops label flipping early without changing the domains", {
    s <- sim_tissue(n = 600)
    full <- hmrf(s$y, s$net, k = 3, beta = 2, zscore = "none", nstart = 10,
        label_tol = 0)
    fast <- hmrf(s$y, s$net, k = 3, beta = 2, zscore = "none",
        init = full$init, blocks = full$blocks)
    expect_lte(fast$iterations, full$iterations)
    expect_gt(ari(fast$class, full$class), 0.95)
    expect_gte(ari(fast$class, s$dom), ari(full$class, s$dom) - 0.01)
})
