.validate_pik <- function(pik, tol = 1e-10) {
  pik <- as.numeric(pik)
  if (length(pik) < 2L || any(!is.finite(pik)) || any(pik <= 0) || any(pik > 1 + tol)) {
    stop("pik must be a finite numeric vector with entries in (0, 1].")
  }
  n <- sum(pik)
  if (abs(n - round(n)) > tol) {
    stop("sum(pik) must be an integer for fixed-size systematic sampling.")
  }
  pik
}

.validate_X <- function(X, N) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != N || ncol(X) < 1L || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix with length(pik) rows and at least one column.")
  }
  X
}

.validate_permutation <- function(sigma, N) {
  sigma <- as.integer(sigma)
  if (length(sigma) != N || !identical(sort(sigma), seq_len(N))) {
    stop("A starting permutation must contain each integer from 1 to N exactly once.")
  }
  sigma
}

.rotate_to_anchor <- function(sigma, anchor = 1L) {
  sigma <- as.integer(sigma)
  k <- match(anchor, sigma)
  if (is.na(k)) stop("anchor is not present in sigma.")
  if (k == 1L) return(sigma)
  c(sigma[k:length(sigma)], sigma[seq_len(k - 1L)])
}

.random_permutation <- function(N) {
  if (N == 2L) return(1:2)
  c(1L, sample.int(N - 1L) + 1L)
}

.balance_scale <- function(X, objective_control) {
  supplied <- objective_control$scale %||% NULL
  if (!is.null(supplied)) {
    supplied <- as.numeric(supplied)
    if (length(supplied) != ncol(X) || any(!is.finite(supplied)) || any(supplied <= 0)) {
      stop("objective_control$scale must contain one positive finite value per column of X.")
    }
    return(supplied)
  }
  ans <- abs(colSums(X))
  if (any(ans < 1e-14)) {
    stop(
      "At least one column total of X is zero (or numerically zero). ",
      "For objective='balance' with metric='relative', supply positive scales through objective_control$scale."
    )
  }
  ans
}

# Prepare the balance criterion.
#
# metric = "relative" reproduces the original criterion
#   sum_j ((Xhat_j - X_j) / a_j)^2
# with a_j = |X_j| by default.
#
# metric = "mahalanobis" implements Cost_2 from Tille (2006, p. 156):
#   (s-pi)' A' (A A')^+ A (s-pi),
# where A = (x_1/pi_1, ..., x_N/pi_N).  Since
# A(s-pi) = Xhat_HT(s) - X, this is a quadratic Mahalanobis/projection
# distance for the imbalance vector.  The book writes an ordinary inverse;
# the Moore-Penrose extension used here coincides with it at full row rank
# and remains well-defined if redundant balancing variables are supplied.
#
# To reuse the highly optimized balance C++ engine without duplicating it,
# construct L such that L L' = (A A')^+ and transform X to X %*% L.
# The existing Euclidean engine with unit scales then evaluates exactly the
# desired quadratic form.
.prepare_balance_control <- function(pik, X, objective_control = list()) {
  metric <- tolower(as.character(objective_control$metric %||% "relative"))[1L]
  metric <- switch(
    metric,
    cost1 = "relative",
    scaled = "relative",
    standardized = "relative",
    cost2 = "mahalanobis",
    projection = "mahalanobis",
    mahal = "mahalanobis",
    metric
  )
  if (!metric %in% c("relative", "mahalanobis")) {
    stop(
      "For objective='balance', objective_control$metric must be ",
      "'relative' or 'mahalanobis'."
    )
  }

  if (metric == "relative") {
    return(list(
      kind = "balance",
      metric = "relative",
      X_engine = X,
      scale = .balance_scale(X, objective_control),
      original_dimension = ncol(X),
      mahalanobis_rank = NA_integer_,
      mahalanobis_tol = NA_real_
    ))
  }

  if (!is.null(objective_control$scale)) {
    warning("objective_control$scale is ignored when metric='mahalanobis'.")
  }
  tol <- as.numeric(objective_control$mahalanobis_tol %||% 1e-10)[1L]
  if (!is.finite(tol) || tol < 0) {
    stop("objective_control$mahalanobis_tol must be one non-negative finite number.")
  }

  # Rows of Z are x_k / pi_k, hence crossprod(Z) = A A'.
  Z <- sweep(X, 1L, pik, "/")
  G <- crossprod(Z)
  eg <- eigen(G, symmetric = TRUE)
  lambda_max <- max(eg$values)
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("The Mahalanobis balance matrix A A' has no positive eigenvalues.")
  }

  cutoff <- tol * lambda_max
  keep <- which(eg$values > cutoff)
  if (!length(keep)) {
    stop("The Mahalanobis balance matrix is numerically rank zero at the requested tolerance.")
  }

  # If G = V diag(lambda) V', then L = V diag(lambda^{-1/2})
  # satisfies L L' = G^+.  X %*% L therefore converts the Mahalanobis
  # criterion into the ordinary squared-Euclidean balance criterion.
  L <- sweep(
    eg$vectors[, keep, drop = FALSE],
    2L,
    sqrt(eg$values[keep]),
    "/"
  )
  X_engine <- X %*% L
  storage.mode(X_engine) <- "double"

  list(
    kind = "balance",
    metric = "mahalanobis",
    X_engine = X_engine,
    scale = rep(1, ncol(X_engine)),
    original_dimension = ncol(X),
    mahalanobis_rank = length(keep),
    mahalanobis_tol = tol,
    mahalanobis_cutoff = cutoff,
    # Keeping the transform is useful for exact diagnostics/validation and
    # costs only p * rank(p) numbers, never an N x N matrix.
    mahalanobis_transform = L
  )
}

# Scale a positive local pik vector to a valid integer expected sample size.
# The truncation min(1, lambda * pik) is found by bisection.  This is used only
# inside the stratified warm start; the final global optimization always uses
# the original pik exactly.
.calibrate_local_pik <- function(pik, target_n, tol = 1e-12) {
  N <- length(pik)
  target_n <- as.integer(target_n)
  if (target_n <= 0L || target_n >= N) return(NULL)

  f <- function(lambda) sum(pmin(1, lambda * pik))
  lo <- 0
  hi <- 1
  while (f(hi) < target_n) hi <- 2 * hi
  for (it in 1:100) {
    mid <- (lo + hi) / 2
    if (f(mid) < target_n) lo <- mid else hi <- mid
    if (hi - lo < tol * max(1, hi)) break
  }
  out <- pmin(1, ((lo + hi) / 2) * pik)
  # Small final correction for floating point arithmetic.
  if (abs(sum(out) - target_n) > 1e-8) {
    stop("Could not calibrate local inclusion probabilities for a stratum.")
  }
  out
}

# ----------------------------------------------------------------------------
# Stratified warm start
# ----------------------------------------------------------------------------

.gcd_int_R <- function(a, b) {
  a <- abs(as.integer(a)); b <- abs(as.integer(b))
  while (b != 0L) { tmp <- a %% b; a <- b; b <- tmp }
  a
}

.make_strata <- function(N, n_strata = NULL, strata = NULL, pik = NULL, tol_pik = 1e-12) {
  if (!is.null(strata)) {
    if (length(strata) != N || anyNA(strata)) stop("strata must have length N and contain no missing values.")
    lev <- unique(strata)
    return(lapply(lev, function(z) which(strata == z)))
  }
  if (is.null(n_strata) || length(n_strata) != 1L || n_strata < 2L || n_strata > N) {
    stop("When stratification=TRUE, supply n_strata between 2 and N, or provide strata explicitly.")
  }
  H <- as.integer(n_strata)

  # For equal pik, construct automatic strata whose sizes are multiples of the
  # reduced systematic period D=N/gcd(N,n).  Then every local stratum has an
  # integer sum of the *original* pik and exactly the same sampling fraction.
  # This aligns the warm-start subproblems with the global systematic phases.
  if (!is.null(pik) && vns_all_equal_pik_cpp(pik, tol_pik)) {
    n <- as.integer(round(sum(pik)))
    D <- N %/% .gcd_int_R(N, n)
    n_blocks <- N %/% D
    if (H <= n_blocks) {
      units <- sample.int(N)
      blocks_per_stratum <- rep(n_blocks %/% H, H)
      if (n_blocks %% H) blocks_per_stratum[seq_len(n_blocks %% H)] <-
        blocks_per_stratum[seq_len(n_blocks %% H)] + 1L
      sizes <- blocks_per_stratum * D
      ends <- cumsum(sizes)
      starts <- c(1L, head(ends, -1L) + 1L)
      return(Map(function(a, b) units[a:b], starts, ends))
    }
  }

  labels <- sample(rep(seq_len(H), length.out = N))
  lapply(seq_len(H), function(h) which(labels == h))
}

.prepare_spread_control <- function(pik, X, objective_control = list()) {
  N <- length(pik)
  metric <- tolower(as.character(objective_control$metric %||% "sb"))[1L]
  if (!metric %in% c("sb", "ib")) stop("For objective='spread', objective_control$metric must be 'sb' or 'ib'.")

  coordinates <- objective_control$coordinates %||% X
  coordinates <- .validate_X(coordinates, N)
  dense_distance_limit <- as.integer(objective_control$dense_distance_limit %||% 2000L)
  if (!is.finite(dense_distance_limit) || dense_distance_limit < 0L) dense_distance_limit <- 0L
  state_cell_limit <- as.numeric(objective_control$state_cell_limit %||% 5e6)
  if (!is.finite(state_cell_limit) || state_cell_limit <= 0) state_cell_limit <- Inf

  wi <- integer(0); wj <- integer(0); wx <- numeric(0)
  W_source <- NULL
  if (metric == "ib") {
    W <- objective_control$W %||% NULL
    if (is.null(W)) {
      if (!requireNamespace("WaveSampling", quietly = TRUE)) {
        stop("metric='ib' requires objective_control$W or package 'WaveSampling'.")
      }
      W <- WaveSampling::wpik(as.matrix(coordinates), pik)
      W_source <- "WaveSampling::wpik"
    } else {
      W_source <- "user"
    }

    if (inherits(W, "sparseMatrix")) {
      if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required for sparse W.")
      sm <- Matrix::summary(W)
      keep <- sm$i != sm$j & is.finite(sm$x) & sm$x != 0
      wi <- as.integer(sm$i[keep]); wj <- as.integer(sm$j[keep]); wx <- as.numeric(sm$x[keep])
    } else {
      W <- as.matrix(W)
      if (!all(dim(W) == c(N, N))) stop("objective_control$W must be an N x N matrix.")
      if (any(!is.finite(W))) stop("W must contain only finite values.")
      diag(W) <- 0
      ij <- which(W != 0, arr.ind = TRUE)
      if (!nrow(ij)) stop("W has no non-zero off-diagonal entries.")
      wi <- as.integer(ij[,1]); wj <- as.integer(ij[,2]); wx <- as.numeric(W[ij])
    }
  }

  list(
    metric = metric,
    coordinates = coordinates,
    wi = wi, wj = wj, wx = wx,
    dense_distance_limit = dense_distance_limit,
    state_cell_limit = state_cell_limit,
    W_source = W_source
  )
}

.objective_setup <- function(pik, X, objective, objective_control = list()) {
  if (objective == "balance") {
    return(.prepare_balance_control(pik, X, objective_control))
  }
  if (objective == "spread") {
    z <- .prepare_spread_control(pik, X, objective_control)
    z$kind <- "spread"
    return(z)
  }
  stop("Unknown objective.")
}

.objective_value_from_setup <- function(sigma, pik, X, objective, setup, tol_pik) {
  sigma <- .validate_permutation(sigma, length(pik))
  if (objective == "balance") {
    return(vns_balance_objective_cpp(
      sigma, pik, setup$X_engine, setup$scale, tol_pik
    ))
  }
  vns_spread_objective_cpp(
    sigma = sigma,
    pik = pik,
    coordinates = setup$coordinates,
    metric = setup$metric,
    wi = setup$wi,
    wj = setup$wj,
    wx = setup$wx,
    tol_pik = tol_pik,
    dense_distance_limit = setup$dense_distance_limit
  )
}

.run_objective_from_setup <- function(sigma, pik, X, objective, setup, search,
                                      improvement_tol, tol_pik, max_iter, max_time) {
  if (objective == "balance") {
    return(vns_run_balance_cpp(
      sigma = sigma, pik = pik, X = setup$X_engine, search = search,
      balance_scale = setup$scale,
      improvement_tol = improvement_tol, tol_pik = tol_pik,
      max_iter = max_iter, max_time = max_time
    ))
  }
  vns_run_spread_cpp(
    sigma = sigma,
    pik = pik,
    coordinates = setup$coordinates,
    metric = setup$metric,
    wi = setup$wi,
    wj = setup$wj,
    wx = setup$wx,
    search = search,
    improvement_tol = improvement_tol,
    tol_pik = tol_pik,
    max_iter = max_iter,
    max_time = max_time,
    dense_distance_limit = setup$dense_distance_limit,
    max_state_cells = setup$state_cell_limit
  )
}

# ----------------------------------------------------------------------------
# Stratified warm start (objective-agnostic high-level logic)
# ----------------------------------------------------------------------------

.stratified_warm_start <- function(pik, X, n_strata, strata,
                                   objective, objective_control, global_setup,
                                   improvement_tol, tol_pik,
                                   strata_control = list(), verbose = FALSE) {
  N <- length(pik)
  groups <- .make_strata(N, n_strata = n_strata, strata = strata, pik = pik, tol_pik = tol_pik)
  H <- length(groups)

  local_max_iter <- as.integer(strata_control$max_iter %||% 500L)
  local_max_time <- as.numeric(strata_control$max_time %||% 1)
  if (local_max_iter < 0L || local_max_time <= 0) stop("Invalid strata_control limits.")

  pieces <- vector("list", H)
  local_info <- vector("list", H)

  for (h in seq_len(H)) {
    units <- groups[[h]]
    Nh <- length(units)
    if (Nh <= 2L) {
      pieces[[h]] <- units[sample.int(Nh)]
      local_info[[h]] <- list(stratum = h, size = Nh, optimized = FALSE)
      next
    }

    ph <- pik[units]
    nh <- as.integer(round(sum(ph)))
    if (nh <= 0L || nh >= Nh) {
      pieces[[h]] <- units[sample.int(Nh)]
      local_info[[h]] <- list(stratum = h, size = Nh, optimized = FALSE)
      next
    }
    ph_work <- if (abs(sum(ph) - nh) <= 1e-10) ph else .calibrate_local_pik(ph, nh)
    Xh <- X[units, , drop = FALSE]

    if (objective == "balance") {
      local_control <- objective_control
      if (identical(global_setup$metric, "relative")) {
        # Preserve the original relative scaling convention in the local
        # warm-start problems.
        local_control$scale <- global_setup$scale * (Nh / N)
      } else {
        # Recompute the Cost_2/Mahalanobis geometry on each local subproblem.
        # This is used only to construct a warm start; the final VNS always
        # optimizes the exact global Mahalanobis criterion.
        local_control$scale <- NULL
        local_control$metric <- "mahalanobis"
        local_control$mahalanobis_tol <- global_setup$mahalanobis_tol
      }
    } else {
      local_control <- objective_control
      local_control$coordinates <- (global_setup$coordinates)[units, , drop = FALSE]
      if (identical(global_setup$metric, "ib")) local_control$W <- NULL
    }
    local_setup <- .objective_setup(ph_work, Xh, objective, local_control)
    sigma0 <- .random_permutation(Nh)
    ans <- .run_objective_from_setup(
      sigma0, ph_work, Xh, objective, local_setup,
      search = "type1",
      improvement_tol = improvement_tol, tol_pik = tol_pik,
      max_iter = local_max_iter, max_time = local_max_time
    )
    pieces[[h]] <- units[ans$sigma]
    local_info[[h]] <- list(
      stratum = h, size = Nh, target_n = nh, optimized = TRUE,
      objective = ans$f, type1_certified = isTRUE(ans$type1_certified),
      stop_reason = ans$stop_reason, elapsed_time = ans$elapsed_time
    )
  }

  if (is.null(strata) && H > 1L) pieces <- pieces[sample.int(H)]
  sigma <- .rotate_to_anchor(unlist(pieces, use.names = FALSE), 1L)
  if (verbose) cat(sprintf("  Stratified warm start built from %d strata.\n", H))
  list(sigma = sigma, strata_info = local_info)
}

# ----------------------------------------------------------------------------
# Public objective helper
# ----------------------------------------------------------------------------

vns_objective <- function(sigma, pik, X, objective = c("balance", "spread"),
                          objective_control = list(), tol_pik = 1e-12) {
  objective <- match.arg(objective)
  pik <- .validate_pik(pik)
  X <- .validate_X(X, length(pik))
  setup <- .objective_setup(pik, X, objective, objective_control)
  .objective_value_from_setup(sigma, pik, X, objective, setup, tol_pik)
}

# TSP is kept as an initialization helper, not as part of the VNS objective.
tsp_order <- function(D, n_starts = 20L, method = "two_opt", seed = NULL,
                      verbose = TRUE) {
  if (!requireNamespace("TSP", quietly = TRUE)) stop("Package 'TSP' is required for tsp_order().")
  D <- as.matrix(D)
  if (nrow(D) != ncol(D) || nrow(D) < 2L || any(!is.finite(D))) stop("D must be a finite square distance matrix.")
  if (!is.null(seed)) set.seed(seed)
  tsp <- TSP::TSP(D)
  best_order <- NULL
  best_length <- Inf
  for (r in seq_len(as.integer(n_starts))) {
    tour <- TSP::solve_TSP(tsp, method = method)
    ord <- as.integer(tour)
    len <- sum(D[cbind(ord, c(ord[-1L], ord[1L]))])
    if (len < best_length) { best_length <- len; best_order <- ord }
  }
  if (verbose) cat(sprintf("TSP warm start: best tour length %.8g (%d starts, method=%s)\n",
                           best_length, as.integer(n_starts), method))
  structure(best_order, tsp_length = best_length, tsp_method = method,
            tsp_n_starts = as.integer(n_starts))
}

# ----------------------------------------------------------------------------
# Main VNS interface
# ----------------------------------------------------------------------------

vns_systematic <- function(
    pik,
    X,
    objective = c("balance", "spread"),
    objective_control = list(),
    start = NULL,
    search = c("adaptive", "type1", "type2"),
    stratification = FALSE,
    n_strata = NULL,
    strata = NULL,
    strata_control = list(),
    n_starts = 1L,
    mixture_top_k = 1L,
    improvement_tol = 1e-10,
    tol_pik = 1e-12,
    max_iter = 100000L,
    max_time = Inf,
    seed = NULL,
    verbose = TRUE
) {
  search <- match.arg(search)
  objective <- match.arg(objective)
  pik <- .validate_pik(pik)
  N <- length(pik)
  X <- .validate_X(X, N)
  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)
  max_time <- as.numeric(max_time)
  mixture_top_k <- as.integer(mixture_top_k)

  if (n_starts < 1L) stop("n_starts must be at least 1.")
  if (mixture_top_k < 1L) stop("mixture_top_k must be at least 1.")
  if (max_iter < 0L) stop("max_iter must be non-negative.")
  if (length(max_time) != 1L || (!is.infinite(max_time) && (!is.finite(max_time) || max_time <= 0))) {
    stop("max_time must be a positive number or Inf.")
  }
  if (improvement_tol < 0 || tol_pik <= 0) stop("Tolerances must be non-negative (tol_pik strictly positive).")
  if (!is.null(seed)) set.seed(seed)

  setup <- .objective_setup(pik, X, objective, objective_control)
  equal_pik <- vns_all_equal_pik_cpp(pik, tol_pik)
  is_mixture <- inherits(start, "systematic_mixture")
  if (isTRUE(stratification) && !is.null(start)) {
    stop("stratification is an initialization strategy and cannot be combined with an explicit start object.")
  }

  strata_metadata <- NULL
  initialization <- NULL
  start_mixture_objective <- NULL
  mixture_component_objectives <- NULL
  source_component <- NULL
  source_weight <- NULL

  if (is_mixture) {
    perms <- start$permutations
    weights0 <- as.numeric(start$weights)
    if (!is.matrix(perms)) perms <- do.call(rbind, perms)
    if (ncol(perms) != N || nrow(perms) != length(weights0)) {
      stop("Invalid systematic_mixture object: dimensions do not match this population.")
    }
    weights0[weights0 < 0 & weights0 > -1e-12] <- 0
    if (any(weights0 < 0) || sum(weights0) <= 0) stop("Mixture weights must be non-negative and have positive sum.")
    weights0 <- weights0 / sum(weights0)
    perms_list <- lapply(seq_len(nrow(perms)), function(i) {
      .rotate_to_anchor(.validate_permutation(perms[i, ], N), 1L)
    })

    # Evaluate mixture components explicitly instead of relying on vapply() here.
    # Besides being easier to debug, this avoids any accidental masking of vapply
    # or argument-matching surprises in user workspaces after sourcing other packages.
    mixture_component_objectives <- numeric(length(perms_list))
    for (i in seq_along(perms_list)) {
      mixture_component_objectives[i] <- .objective_value_from_setup(
        sigma = perms_list[[i]],
        pik = pik,
        X = X,
        objective = objective,
        setup = setup,
        tol_pik = tol_pik
      )
    }
    start_mixture_objective <- sum(weights0 * mixture_component_objectives)

    # Expected sample cost is linear in the sampling design. Thus a convex
    # mixture cannot improve on its best component for the chosen objective.
    ord <- order(mixture_component_objectives)
    keep <- head(ord, min(mixture_top_k, length(ord)))
    starts <- perms_list[keep]
    source_component <- keep
    source_weight <- weights0[keep]
    initialization <- if (length(keep) == 1L) "systematic_mixture_best_component" else "systematic_mixture_multistart"
    if (n_starts != 1L) warning("n_starts is ignored for a systematic_mixture start; use mixture_top_k instead.")
  } else if (!is.null(start)) {
    starts <- list(.rotate_to_anchor(.validate_permutation(start, N), 1L))
    initialization <- "user_permutation"
    if (n_starts != 1L) warning("n_starts is ignored when an explicit starting permutation is supplied.")
  } else {
    starts <- vector("list", n_starts)
    strata_metadata <- vector("list", n_starts)
    for (r in seq_len(n_starts)) {
      if (isTRUE(stratification)) {
        sws <- .stratified_warm_start(
          pik, X, n_strata, strata,
          objective = objective, objective_control = objective_control,
          global_setup = setup,
          improvement_tol = improvement_tol, tol_pik = tol_pik,
          strata_control = strata_control,
          verbose = verbose && r == 1L
        )
        starts[[r]] <- sws$sigma
        strata_metadata[[r]] <- sws$strata_info
      } else {
        starts[[r]] <- .random_permutation(N)
      }
    }
    initialization <- if (isTRUE(stratification)) "stratified_warm_start" else "random"
  }

  if (verbose) {
    cat("VNS systematic sampling optimization\n")
    cat("------------------------------------\n")
    cat(sprintf("N: %d | n: %d | objective: %s", N, round(sum(pik)), objective))
    if (objective %in% c("balance", "spread")) cat(sprintf(" (%s)", toupper(setup$metric)))
    cat("\n")
    cat(sprintf("pik: %s | engine: C++ %s\n",
                if (equal_pik) "equal" else "unequal",
                if (equal_pik) "incremental" else "general"))
    cat(sprintf("initialization: %s | search: %s\n", initialization, search))
    if (is_mixture) {
      cat(sprintf("mixture components: %d | optimized warm starts: %d\n", nrow(perms), length(starts)))
      cat(sprintf("mixture objective: %.8g | best component: %.8g\n",
                  start_mixture_objective, min(mixture_component_objectives)))
    }
  }

  runs <- vector("list", length(starts))
  initial_objectives <- numeric(length(starts))
  for (r in seq_along(starts)) {
    if (verbose && length(starts) > 1L) cat(sprintf("  Optimizing start %d of %d...\n", r, length(starts)))
    ans <- .run_objective_from_setup(
      starts[[r]], pik, X, objective, setup, search,
      improvement_tol, tol_pik, max_iter, max_time
    )
    initial_objectives[r] <- as.numeric(ans$initial_f)
    numerical_guard <- max(1e-9, 100 * improvement_tol)
    if (ans$f > initial_objectives[r] + numerical_guard) {
      stop("Internal monotonicity check failed: the optimized design is worse than its start.")
    }
    runs[[r]] <- ans
  }

  fvals <- numeric(length(runs))
  for (r in seq_along(runs)) fvals[r] <- as.numeric(runs[[r]]$f)
  best <- which.min(fvals)
  final_perms <- matrix(runs[[best]]$sigma, nrow = 1L)
  final_obj <- fvals[best]
  if (is_mixture) initial_obj <- min(mixture_component_objectives) else initial_obj <- initial_objectives[best]
  chosen <- best

  diag_rows <- lapply(seq_along(runs), function(r) {
    z <- runs[[r]]
    rel <- if (abs(initial_objectives[r]) > 1e-12) (initial_objectives[r] - z$f) / abs(initial_objectives[r]) else NA_real_
    data.frame(
      run = r,
      source_component = if (is_mixture) source_component[r] else NA_integer_,
      source_weight = if (is_mixture) source_weight[r] else NA_real_,
      initial_objective = initial_objectives[r], final_objective = z$f,
      improvement = initial_objectives[r] - z$f, relative_improvement = rel,
      iterations = z$iterations, elapsed_time = z$elapsed_time,
      evaluation_method = z$evaluation_method,
      type1_certified = z$type1_certified, type2_certified = z$type2_certified,
      local_minimum_type = z$local_minimum_type, stop_reason = z$stop_reason,
      n_adjacent_evaluations = z$n_adjacent_evaluations,
      n_full_evaluations = z$n_full_evaluations,
      n_adjacent_moves = z$n_adjacent_moves, n_global_moves = z$n_global_moves,
      n_type1_certificates = z$n_type1_certificates,
      stringsAsFactors = FALSE
    )
  })
  diagnostics <- do.call(rbind, diag_rows)

  relative_improvement <- if (abs(initial_obj) > 1e-12) (initial_obj - final_obj) / abs(initial_obj) else NA_real_
  control_out <- if (objective == "balance") {
    if (identical(setup$metric, "relative")) {
      list(metric = "relative", scale = setup$scale)
    } else {
      list(
        metric = "mahalanobis",
        mahalanobis_rank = setup$mahalanobis_rank,
        original_dimension = setup$original_dimension,
        mahalanobis_tol = setup$mahalanobis_tol
      )
    }
  } else {
    list(metric = setup$metric, coordinates = setup$coordinates,
         dense_distance_limit = setup$dense_distance_limit,
         state_cell_limit = setup$state_cell_limit, W_source = setup$W_source)
  }

  out <- list(
    call = match.call(), objective = objective, objective_control = control_out,
    N = N, n = round(sum(pik)), pik = pik, equal_pik = equal_pik,
    evaluation_method = runs[[best]]$evaluation_method,
    initialization = initialization, search = search,
    initial_objective = initial_obj, final_objective = final_obj,
    improvement = initial_obj - final_obj, relative_improvement = relative_improvement,
    permutations = final_perms, weights = 1, sigma = as.integer(final_perms[1, ]),
    diagnostics = diagnostics, runs = runs, selected_runs = chosen,
    strata = strata_metadata,
    start_object = if (is_mixture) start else NULL,
    start_mixture_objective = start_mixture_objective,
    mixture_component_objectives = mixture_component_objectives
  )
  class(out) <- "vns_systematic_fit"
  if (verbose) print(out)
  out
}

print.vns_systematic_fit <- function(x, ...) {
  cat("\nVNS result\n")
  label <- x$objective
  if (identical(x$objective, "spread")) {
    label <- sprintf("spread (%s)", toupper(x$objective_control$metric))
  } else if (identical(x$objective, "balance")) {
    label <- sprintf("balance (%s)", toupper(x$objective_control$metric %||% "relative"))
  }
  cat(sprintf("  objective:          %s\n", label))
  cat(sprintf("  initial objective:  %.8g\n", x$initial_objective))
  cat(sprintf("  final objective:    %.8g\n", x$final_objective))
  cat(sprintf("  improvement:        %.8g\n", x$improvement))
  if (is.finite(x$relative_improvement))
    cat(sprintf("  relative reduction: %.2f%%\n", 100 * x$relative_improvement))
  cat(sprintf("  engine:             %s\n", x$evaluation_method))
  if (!is.null(x$start_mixture_objective))
    cat(sprintf("  original mixture:   %.8g (best component %.8g)\n",
                x$start_mixture_objective, min(x$mixture_component_objectives)))

  d <- x$diagnostics[x$selected_runs, , drop = FALSE]
  cat(sprintf("  local minimum:      %s\n", d$local_minimum_type))
  cat(sprintf("  stop reason:        %s\n", d$stop_reason))
  cat(sprintf("  accepted moves:     %d (%d adjacent, %d global)\n",
              d$iterations, d$n_adjacent_moves, d$n_global_moves))
  cat(sprintf("  elapsed time:       %.3f s\n", d$elapsed_time))
  invisible(x)
}

summary.vns_systematic_fit <- function(object, ...) {
  object$diagnostics
}

plot.vns_systematic_fit <- function(x, component = 1L, ...) {
  component <- as.integer(component)
  if (component < 1L || component > length(x$runs)) stop("Invalid component/run index.")
  path <- x$runs[[component]]$f_path
  positive <- all(is.finite(path)) && all(path > 0)
  plot(seq_along(path) - 1L, path, type = "l", log = if (positive) "y" else "",
       xlab = "Accepted move",
       ylab = if (positive) "Objective (log scale)" else "Objective", ...)
  invisible(x)
}

# ----------------------------------------------------------------------------
# Generic Monte Carlo estimation of a sampling design
# ----------------------------------------------------------------------------

.normalize_sample <- function(s, N) {
  if (is.list(s)) {
    if (!is.null(s$sample)) s <- s$sample
    else if (!is.null(s$s)) s <- s$s
  }
  if ((is.numeric(s) || is.logical(s)) && length(s) == N) {
    sr <- round(as.numeric(s))
    if (all(is.finite(s)) && all(abs(as.numeric(s) - sr) <= 1e-8) && all(sr %in% c(0, 1))) {
      return(as.integer(sr))
    }
  }
  if (is.numeric(s) && length(s) >= 1L && all(is.finite(s)) && all(s == as.integer(s)) &&
      all(s >= 1L) && all(s <= N) && length(unique(s)) == length(s)) {
    out <- integer(N)
    out[as.integer(s)] <- 1L
    return(out)
  }
  stop("sampler() must return either a 0/1 vector of length N or the selected unit indices.")
}

.sample_pair_set <- function(N, K) {
  total <- N * (N - 1) / 2
  K <- min(as.integer(K), total)
  if (K <= 0L) stop("n_pairs must be positive.")
  if (total <= 20000 && K == total) return(t(utils::combn(N, 2L)))

  ans <- matrix(integer(0), ncol = 2L)
  while (nrow(ans) < K) {
    need <- K - nrow(ans)
    batch <- max(2L * need, 1000L)
    a <- sample.int(N, batch, replace = TRUE)
    b <- sample.int(N, batch, replace = TRUE)
    keep <- a != b
    lo <- pmin(a[keep], b[keep])
    hi <- pmax(a[keep], b[keep])
    cand <- unique(cbind(lo, hi))
    ans <- unique(rbind(ans, cand))
  }
  ans[seq_len(K), , drop = FALSE]
}

estimate_sampling_design_mc <- function(sampler, N, M = 2000L,
                                        n_pairs = 2000L, pairs = NULL,
                                        seed = NULL, verbose = TRUE) {
  if (!is.function(sampler)) stop("sampler must be a function.")
  N <- as.integer(N); M <- as.integer(M)
  if (N < 2L || M < 2L) stop("N and M must be at least 2.")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(pairs)) pairs <- .sample_pair_set(N, n_pairs)
  pairs <- as.matrix(pairs)
  storage.mode(pairs) <- "integer"
  if (ncol(pairs) != 2L || any(pairs < 1L) || any(pairs > N) || any(pairs[,1] == pairs[,2])) {
    stop("pairs must be a valid two-column matrix of population indices.")
  }

  pik_sum <- numeric(N)
  pi2_sum <- numeric(nrow(pairs))
  sample_sizes <- integer(M)
  sampler_has_arg <- length(formals(sampler)) > 0L

  if (verbose) cat(sprintf("Estimating design moments from %d Monte Carlo samples...\n", M))
  for (m in seq_len(M)) {
    raw <- if (sampler_has_arg) sampler(m) else sampler()
    s <- .normalize_sample(raw, N)
    sample_sizes[m] <- sum(s)
    pik_sum <- pik_sum + s
    pi2_sum <- pi2_sum + s[pairs[,1]] * s[pairs[,2]]
  }

  out <- list(
    N = N,
    M = M,
    pik_hat = pik_sum / M,
    pairs = pairs,
    pi2_hat = pi2_sum / M,
    sample_size = list(
      mean = mean(sample_sizes),
      sd = stats::sd(sample_sizes),
      min = min(sample_sizes),
      max = max(sample_sizes)
    )
  )
  class(out) <- "empirical_sampling_design"
  out
}

# A CUBE-specific simulator is deliberately only a thin wrapper.  The generic
# estimator above knows nothing about CUBE and can be reused with any sampler.
estimate_cube_design <- function(pik, X_balance, M = 2000L,
                                 n_pairs = 2000L, pairs = NULL,
                                 force_fixed_size = TRUE,
                                 seed = NULL, verbose = TRUE) {
  if (!requireNamespace("sampling", quietly = TRUE)) {
    stop("Package 'sampling' is required for estimate_cube_design().")
  }
  pik <- .validate_pik(pik)
  N <- length(pik)
  X_balance <- .validate_X(X_balance, N)
  # In samplecube(), balancing on x_k = pik_k gives
  # sum_{k in s} x_k/pik_k = |s| = sum(pik), hence an exactly fixed
  # sample size (up to the landing tolerance).  This makes comparisons with
  # fixed-size systematic designs fair by default.
  X_cube <- if (isTRUE(force_fixed_size)) cbind(.fixed_size = pik, X_balance) else X_balance

  eps <- 1e-12
  fixed_in <- pik >= 1 - eps
  fixed_out <- pik <= eps
  active <- !(fixed_in | fixed_out)

  sampler <- function() {
    s <- integer(N)
    s[fixed_in] <- 1L
    if (any(active)) {
      res <- sampling::samplecube(
        X_cube[active, , drop = FALSE],
        pik[active], order = 1, comment = FALSE
      )
      s[active] <- .normalize_sample(res, sum(active))
    }
    s
  }
  out <- estimate_sampling_design_mc(
    sampler = sampler, N = N, M = M, n_pairs = n_pairs,
    pairs = pairs, seed = seed, verbose = verbose
  )
  out$generator <- "cube"
  out$force_fixed_size <- isTRUE(force_fixed_size)
  if (isTRUE(force_fixed_size)) {
    n_target <- round(sum(pik))
    if (out$sample_size$min != n_target || out$sample_size$max != n_target) {
      warning("CUBE Monte Carlo samples were not all of the requested fixed size. Check the balancing matrix and landing behavior.")
    }
  }
  out
}

# ----------------------------------------------------------------------------
# Convex approximation by systematic designs
# ----------------------------------------------------------------------------

.project_simplex <- function(v) {
  if (length(v) == 1L) return(1)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  idx <- which(u - cssv / seq_along(u) > 0)
  if (!length(idx)) return(rep(1 / length(v), length(v)))
  rho <- max(idx)
  theta <- cssv[rho] / rho
  pmax(v - theta, 0)
}

.largest_eigenvalue_power <- function(G, n_iter = 100L) {
  p <- nrow(G)
  if (p == 1L) return(as.numeric(G[1,1]))
  v <- rep(1 / sqrt(p), p)
  for (it in seq_len(n_iter)) {
    w <- as.numeric(G %*% v)
    nw <- sqrt(sum(w * w))
    if (nw <= .Machine$double.eps) return(0)
    v <- w / nw
  }
  as.numeric(crossprod(v, G %*% v))
}

.fit_simplex_least_squares <- function(A, b, max_iter = 5000L, tol = 1e-9) {
  A <- as.matrix(A); b <- as.numeric(b)
  R <- ncol(A)
  G <- crossprod(A)
  cvec <- as.numeric(crossprod(A, b))
  L <- .largest_eigenvalue_power(G)
  if (!is.finite(L) || L <= .Machine$double.eps) {
    return(list(alpha = rep(1 / R, R), iterations = 0L, converged = TRUE))
  }

  alpha <- rep(1 / R, R)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    grad <- as.numeric(G %*% alpha - cvec)
    new_alpha <- .project_simplex(alpha - grad / L)
    if (max(abs(new_alpha - alpha)) <= tol) {
      alpha <- new_alpha
      converged <- TRUE
      break
    }
    alpha <- new_alpha
  }
  list(alpha = alpha, iterations = it, converged = converged)
}

.as_permutation_matrix <- function(x, N) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.matrix(x)) x <- do.call(rbind, x)
  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  if (ncol(x) != N) stop("candidate_permutations must have N columns.")
  for (r in seq_len(nrow(x))) x[r, ] <- .rotate_to_anchor(.validate_permutation(x[r, ], N), 1L)
  x
}

approximate_design_by_systematics <- function(
    design,
    pik,
    n_candidates = 100L,
    candidate_permutations = NULL,
    max_components = 20L,
    alpha_tol = 1e-8,
    regression_max_iter = 5000L,
    regression_tol = 1e-9,
    tol_pik = 1e-12,
    seed = NULL,
    verbose = TRUE
) {
  if (!inherits(design, "empirical_sampling_design")) {
    stop("design must be an object returned by estimate_sampling_design_mc() or estimate_cube_design().")
  }
  pik <- .validate_pik(pik)
  N <- length(pik)
  if (design$N != N) stop("design and pik have incompatible population sizes.")
  target_n <- sum(pik)
  if (!is.null(design$sample_size)) {
    fixed_mc <- identical(design$sample_size$min, design$sample_size$max)
    if (!fixed_mc || abs(design$sample_size$mean - target_n) > 1e-8) {
      warning(
        "The Monte Carlo source design is not fixed-size with sample size sum(pik). ",
        "The systematic mixture will still have the supplied pik and fixed size, so the approximation is only moment-based."
      )
    }
  }
  if (!is.null(seed)) set.seed(seed)

  n_candidates <- as.integer(n_candidates)
  max_components <- as.integer(max_components)
  if (n_candidates < 1L || max_components < 1L) stop("n_candidates and max_components must be positive.")

  supplied <- .as_permutation_matrix(candidate_permutations, N)
  n_random <- max(0L, n_candidates - if (is.null(supplied)) 0L else nrow(supplied))
  random_perms <- if (n_random > 0L) {
    do.call(rbind, replicate(n_random, .random_permutation(N), simplify = FALSE))
  } else NULL
  perms <- if (is.null(supplied)) random_perms else if (is.null(random_perms)) supplied else rbind(supplied, random_perms)
  if (is.null(perms) || nrow(perms) < 1L) stop("No candidate systematic designs were generated.")

  pairs <- design$pairs
  target <- as.numeric(design$pi2_hat)
  K <- nrow(pairs)
  R <- nrow(perms)
  A <- matrix(0, nrow = K, ncol = R)

  if (verbose) {
    cat(sprintf("Approximating the design with %d candidate systematic designs\n", R))
    cat(sprintf("using %d second-order inclusion probabilities...\n", K))
  }
  for (r in seq_len(R)) {
    A[, r] <- systematic_pair_probs_cpp(perms[r, ], pik, pairs, tol_pik)
    if (verbose && (r %% max(1L, floor(R / 10L)) == 0L || r == R)) {
      cat(sprintf("  candidates evaluated: %d/%d\n", r, R))
    }
  }

  fit <- .fit_simplex_least_squares(A, target, regression_max_iter, regression_tol)
  alpha <- fit$alpha

  # Sparsify only after the convex fit, then refit on the retained columns so
  # that the weights still solve a simplex-constrained problem.
  keep <- which(alpha > alpha_tol)
  if (!length(keep)) keep <- which.max(alpha)
  if (length(keep) > max_components) keep <- order(alpha, decreasing = TRUE)[seq_len(max_components)]
  if (length(keep) < R) {
    fit2 <- .fit_simplex_least_squares(A[, keep, drop = FALSE], target,
                                      regression_max_iter, regression_tol)
    alpha2 <- fit2$alpha
  } else {
    alpha2 <- alpha
    fit2 <- fit
  }

  fitted <- as.numeric(A[, keep, drop = FALSE] %*% alpha2)
  rmse <- sqrt(mean((fitted - target)^2))
  max_abs <- max(abs(fitted - target))

  out <- list(
    permutations = perms[keep, , drop = FALSE],
    weights = as.numeric(alpha2 / sum(alpha2)),
    pairs = pairs,
    target_pi2 = target,
    fitted_pi2 = fitted,
    diagnostics = list(
      n_candidates = R,
      n_components = length(keep),
      pairwise_rmse = rmse,
      pairwise_max_abs_error = max_abs,
      regression_iterations = fit2$iterations,
      regression_converged = fit2$converged,
      mc_first_order_max_error = max(abs(design$pik_hat - pik))
    ),
    source_design = design
  )
  class(out) <- "systematic_mixture"
  if (verbose) print(out)
  out
}

print.systematic_mixture <- function(x, ...) {
  cat("\nSystematic-mixture approximation\n")
  cat(sprintf("  components:              %d\n", length(x$weights)))
  cat(sprintf("  pairwise RMSE:            %.6g\n", x$diagnostics$pairwise_rmse))
  cat(sprintf("  max pairwise abs. error:  %.6g\n", x$diagnostics$pairwise_max_abs_error))
  cat(sprintf("  max MC first-order error: %.6g\n", x$diagnostics$mc_first_order_max_error))
  cat("  note: when passed to vns_systematic(), components are treated as warm starts;\n")
  cat("        the final optimized result is a single systematic design.\n")
  invisible(x)
}
