#' @title Spatial domains with a hidden Markov random field
#' @name hmrf
#' @description
#' Fits the hidden Markov random field (HMRF) model of Zhu et al. (2018), as
#' in the \pkg{smfishHmrf} package: each cell's values follow a multivariate
#' Gaussian per domain, and a Potts prior of strength `beta` favours domains
#' shared with spatial neighbours. Fitted by expectation-maximization with
#' iterated conditional modes (ICM) over the colour classes of the network,
#' from a k-means start. Each colour class is an independent set, so its ICM
#' update and the neighbour counts are sparse matrix products; the results
#' equal \pkg{smfishHmrf}'s on the same inputs.
#'
#' **If you use it, please cite the original paper** (see References).
#' @param y cells x features numeric matrix, cell IDs as row names (e.g. the
#'   scaled expression of spatial genes)
#' @param network spatial network: a `data.frame` with columns `from` and `to`
#'   (cell IDs), each pair listed once
#' @param k integer, number of domains
#' @param beta numeric, strength of the spatial prior (0: no prior)
#' @param zscore `"cell"` (default, as Giotto's `initHMRF_V2()`): z-score
#'   each cell across the features; `"none"`: use `y` as given
#' @param init optional start from [hmrf_init()], e.g. to reuse one start for
#'   several `beta`; computed from `y` when `NULL`
#' @param blocks optional colour classes from [hmrf_blocks()]; computed from
#'   `network` when `NULL`
#' @param nstart,seed k-means starts and seed (used when `init` is `NULL`)
#' @param maxit integer, maximum EM iterations
#' @param err numeric, relative change of means and covariances at which EM
#'   stops
#' @param label_tol numeric, share of cells changing domain at which EM stops
#'   (labels can flip between two states without the parameters converging);
#'   0 stops on `err` or `maxit` only, as \pkg{smfishHmrf}
#' @param n_icm integer, ICM sweeps per EM iteration
#' @param damp_tolerance,damp_factor see [hmrf_init()]
#' @returns A list: `class` (integer domain per cell), `prob` (cells x k
#'   posterior probabilities), `mu` (features x k), `sigma` (features x
#'   features x k), `likelihood`, `iterations`, `init` and `blocks`.
#' @references Zhu Q, Shah S, Dries R, Cai L, Yuan GC. Identification of
#'   spatially associated subpopulations by combining scRNAseq and sequential
#'   fluorescence in situ hybridization data. *Nature Biotechnology* 36,
#'   1183–1190 (2018). \doi{10.1038/nbt.4260}
#' @examples
#' set.seed(1)
#' xy <- matrix(runif(400), ncol = 2) * 10
#' dom <- 1L + (xy[, 1] > 5)
#' y <- matrix(rnorm(200 * 4), 200) + c(-1.5, 1.5)[dom]
#' rownames(y) <- paste0("c", 1:200)
#' nn <- function(xy, k) {
#'     d <- as.matrix(stats::dist(xy)); diag(d) <- Inf
#'     t(apply(d, 1, order))[, seq_len(k)]
#' }
#' knn <- nn(xy, 5)
#' e <- unique(t(apply(cbind(rep(1:200, 5), as.vector(knn)), 1, sort)))
#' net <- data.frame(from = rownames(y)[e[, 1]], to = rownames(y)[e[, 2]])
#' fit <- hmrf(y, net, k = 2, beta = 2, zscore = "none", nstart = 10)
#' table(fit$class, dom)
#' @md
#' @export
hmrf <- function(y, network, k, beta = 2, zscore = c("cell", "none"),
    init = NULL, blocks = NULL, nstart = 1000L, seed = 100L, maxit = 50L,
    err = 1e-7, label_tol = 1e-4, n_icm = 10L, damp_tolerance = 1e-5,
    damp_factor = 1.05) {
    zscore <- match.arg(zscore)
    if (is.null(rownames(y))) {
        stop("[hmrf] `y` needs cell IDs as row names", call. = FALSE)
    }
    stopifnot(k >= 2, beta >= 0, maxit >= 1, n_icm >= 1, label_tol >= 0)
    y <- as.matrix(y)
    if (identical(zscore, "cell")) y <- t(scale(t(y)))
    y[!is.finite(y)] <- 0
    adj <- hmrf_adjacency(network, rownames(y))
    if (is.null(blocks)) blocks <- hmrf_blocks(adj)
    if (is.null(init)) {
        init <- hmrf_init(y, k, nstart = nstart, seed = seed,
            damp_tolerance = damp_tolerance, damp_factor = damp_factor)
    }
    fit <- .hmrf_fit(y, adj, blocks, init, k = as.integer(k), beta = beta,
        maxit = maxit, err = err, label_tol = label_tol, n_icm = n_icm)
    rownames(fit$prob) <- rownames(y)
    c(fit, list(init = init, blocks = blocks))
}

#' @title Network helpers
#' @name hmrf_adjacency
#' @description `hmrf_adjacency()` turns an edge list of cell IDs into a
#'   symmetric sparse matrix of neighbour counts (self edges dropped);
#'   `hmrf_blocks()` colours it with igraph's DSatur heuristic and returns the
#'   colour classes, which are independent sets.
#' @param network `data.frame` with columns `from` and `to` (cell IDs)
#' @param ids cell IDs, the row order of `y`
#' @param adj a matrix from `hmrf_adjacency()`
#' @returns `hmrf_adjacency()`: a sparse `n x n` matrix. `hmrf_blocks()`: a
#'   list of integer vectors (row positions) partitioning the cells.
#' @export
hmrf_adjacency <- function(network, ids) {
    i <- match(as.character(network$from), ids)
    j <- match(as.character(network$to), ids)
    keep <- !is.na(i) & !is.na(j) & i != j
    i <- i[keep]
    j <- j[keep]
    n <- length(ids)
    Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(n, n))
}

#' @rdname hmrf_adjacency
#' @export
hmrf_blocks <- function(adj) {
    g <- igraph::graph_from_adjacency_matrix(adj > 0, mode = "undirected",
        diag = FALSE)
    col <- igraph::greedy_vertex_coloring(g, heuristic = "dsatur")
    lapply(sort(unique(col)), function(cl) which(col == cl))
}

#' @title HMRF start
#' @name hmrf_init
#' @description k-means start, per-domain covariances and damping, as in
#'   \pkg{smfishHmrf}: a domain whose covariance determinant is below
#'   `damp_tolerance` gets the smallest diagonal term (searched in steps of
#'   `damp_factor` from 1e-4) that brings it above.
#' @param y cells x features matrix (already z-scored if wanted)
#' @param k integer, number of domains
#' @param nstart integer, k-means starts
#' @param seed integer, seed (0 or negative leaves the RNG alone)
#' @param damp_tolerance,damp_factor damping search settings
#' @returns A list with `mu` (features x k), `sigma` (features x features x
#'   k) and `damp` (length k).
#' @export
hmrf_init <- function(y, k, nstart = 1000L, seed = 100L,
    damp_tolerance = 1e-5, damp_factor = 1.05) {
    if (seed > 0) set.seed(seed)
    kk <- stats::kmeans(y, k, nstart = nstart, iter.max = 100)
    m <- ncol(y)
    sigma <- array(0, c(m, m, k))
    damp <- numeric(k)
    for (j in seq_len(k)) {
        sigma[, , j] <- stats::cov(y[kk$cluster == j, , drop = FALSE])
        d <- .hmrf_damp(sigma[, , j], damp_factor, damp_tolerance)
        damp[j] <- if (is.null(d)) 0 else d
    }
    list(mu = t(kk$centers), sigma = sigma, damp = damp)
}

.hmrf_damp <- function(sigma, factor, cutoff, start = 1e-4) {
    if (det(sigma) >= cutoff) return(NULL)
    eye <- diag(nrow(sigma))
    damp <- start
    if (det(sigma + eye * damp) < cutoff) {
        repeat {
            damp <- damp * factor
            if (det(sigma + eye * damp) > cutoff) break
        }
    } else {
        repeat {
            damp <- damp / factor
            if (det(sigma + eye * damp) < cutoff) {
                damp <- damp * factor
                break
            }
        }
    }
    damp
}

# Gaussian densities of every cell under every domain (clipped as smfishHmrf)
.hmrf_density <- function(y, mu, sigma, damp) {
    n <- nrow(y)
    m <- ncol(y)
    k <- ncol(mu)
    den <- matrix(0, n, k)
    for (j in seq_len(k)) {
        r <- chol(sigma[, , j] + diag(m) * damp[j])
        ri <- backsolve(r, diag(m))
        yr <- y %*% ri
        c0 <- drop(mu[, j] %*% ri)
        maha <- rowSums(yr * yr) - 2 * drop(yr %*% c0) + sum(c0 * c0)
        e <- -0.5 * maha - ((m / 2) * log(2 * pi) + sum(log(diag(r))))
        den[, j] <- exp(pmin(700, pmax(e, -700)))
    }
    den
}

.hmrf_onehot <- function(cl, k) {
    z <- matrix(0, length(cl), k)
    z[cbind(seq_along(cl), cl)] <- 1
    z
}

.hmrf_relerr <- function(a, b) max(abs(a - b)) / (1 + max(abs(a), abs(b)))

.hmrf_fit <- function(y, adj, blocks, init, k, beta, maxit, err, n_icm,
    label_tol = 0) {
    mu <- init$mu
    sigma <- init$sigma
    damp <- init$damp
    z <- .hmrf_onehot(max.col(.hmrf_density(y, mu, sigma, damp)), k)
    it <- 0L
    cl_old <- max.col(z, "first")
    repeat {
        it <- it + 1L
        mu_old <- mu
        sigma_old <- sigma
        den <- .hmrf_density(y, mu, sigma, damp)
        logden <- log(den)
        for (s in seq_len(n_icm)) {
            for (b in blocks) {
                ni <- as.matrix(adj[b, , drop = FALSE] %*% z)
                cl <- max.col(logden[b, , drop = FALSE] + beta * ni,
                    ties.method = "first")
                z[b, ] <- .hmrf_onehot(cl, k)
            }
        }
        nij <- as.matrix(adj %*% z)
        prob <- den * exp(beta * nij)
        prob <- prob / rowSums(prob)
        prob[!is.finite(rowSums(prob)), ] <- 1 / k
        unnorm <- den * exp(pmin(700, pmax(beta * nij, -700)))
        w <- colSums(prob)
        mu <- crossprod(y, prob) %*% diag(1 / w, k)
        for (j in seq_len(k)) {
            sigma[, , j] <- crossprod(y * sqrt(prob[, j])) / w[j] -
                tcrossprod(mu[, j])
        }
        cl <- max.col(prob, "first")
        flipped <- mean(cl != cl_old)
        cl_old <- cl
        if ((.hmrf_relerr(mu_old, mu) < err &&
            .hmrf_relerr(sigma_old, sigma) < err) || flipped < label_tol ||
            it > maxit) break
    }
    t_prob <- colMeans(unnorm / rowSums(unnorm))
    list(class = max.col(prob, "first"), prob = prob, mu = mu, sigma = sigma,
        likelihood = sum(log(unnorm %*% t_prob)), iterations = it)
}
