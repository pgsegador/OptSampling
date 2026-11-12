.cos_validate_problem <- function(values, pik) {
  values <- as.numeric(values)
  pik <- as.numeric(pik)

  if (length(values) != length(pik)) {
    stop("values and pik must have the same length.")
  }
  if (length(values) < 2L) stop("At least two population units are required.")
  if (any(!is.finite(values))) stop("values must be finite.")
  if (any(!is.finite(pik)) || any(pik <= 0) || any(pik >= 1)) {
    stop("All inclusion probabilities must lie strictly between 0 and 1.")
  }

  n <- as.integer(round(sum(pik)))
  if (abs(sum(pik) - n) > 1e-8) {
    stop("sum(pik) must be an integer for the fixed-size problem.")
  }
  if (n <= 0L || n >= length(pik)) stop("The sample size must satisfy 0 < n < N.")

  list(values = values, pik = pik, N = length(pik), n = n)
}


.cos_sample_key <- function(s) paste(sort(as.integer(s)), collapse = ",")


.cos_basis_matrix <- function(samples, N) {
  samples <- as.matrix(samples)
  storage.mode(samples) <- "integer"
  m <- nrow(samples)
  n <- ncol(samples)

  Matrix::sparseMatrix(
    i = as.integer(as.vector(t(samples))),
    j = rep(seq_len(m), each = n),
    x = 1,
    dims = c(N, m),
    giveCsparse = TRUE
  )
}


.cos_sample_sums <- function(samples, z) {
  samples <- as.matrix(samples)
  m <- nrow(samples)
  zz <- matrix(
    z[as.integer(as.vector(t(samples)))],
    nrow = m,
    byrow = TRUE
  )
  rowSums(zz)
}


.cos_sample_costs <- function(samples, z, target) {
  (.cos_sample_sums(samples, z) - target)^2
}


.cos_choose_within <- function(N, n, threshold) {
  lchoose(N, n) <= log(threshold)
}


# -----------------------------------------------------------------------------
# Exact-feasible initial basis
# -----------------------------------------------------------------------------
# A systematic design provides a convex decomposition of pik into fixed-size
# sample indicators.  For equal or otherwise non-generic pik the systematic
# support can contain fewer than N samples, which is a degenerate BFS.
#
# We obtain a full N-column basis by applying an arbitrarily small zero-sum
# perturbation ONLY to discover a neighbouring systematic simplex.  Crucially,
# the probabilities of the basic variables are then recomputed from
#
#        B x = pik
#
# using the ORIGINAL pik.  Thus the prescribed inclusion probabilities are
# never perturbed.  Zero basic variables naturally appear in degenerate cases.

.cos_initial_basis <- function(
    pik,
    n,
    perturb_scale = 1e-7,
    max_tries = 30L,
    feasibility_tol = 1e-9,
    seed = NULL,
    verbose = FALSE) {

  N <- length(pik)
  if (!is.null(seed)) set.seed(seed)

  eps_bound <- 1e-12

  for (attempt in 0:max_tries) {
    if (attempt == 0L) {
      q <- pik
      sigma <- seq_len(N)
      step <- 0
    } else {
      direction <- rnorm(N)
      direction <- direction - mean(direction)
      md <- max(abs(direction))
      if (!is.finite(md) || md == 0) next
      direction <- direction / md

      bounds <- numeric(0)
      neg <- direction < 0
      pos <- direction > 0
      if (any(neg)) {
        bounds <- c(bounds, (pik[neg] - eps_bound) / (-direction[neg]))
      }
      if (any(pos)) {
        bounds <- c(bounds, ((1 - eps_bound) - pik[pos]) / direction[pos])
      }
      max_step <- min(bounds)
      if (!is.finite(max_step) || max_step <= 0) next

      step <- min(perturb_scale, 0.05 * max_step)
      q <- pik + step * direction
      # Correct only floating-point drift in the zero-sum perturbation.
      q <- q + (sum(pik) - sum(q)) / N
      sigma <- sample.int(N)
    }

    supp <- tryCatch(
      cos_systematic_support_cpp(q, as.integer(sigma), tol = 1e-13),
      error = function(e) NULL
    )
    if (is.null(supp)) next

    samples <- as.matrix(supp$samples)
    if (nrow(samples) != N || ncol(samples) != n) next

    B <- .cos_basis_matrix(samples, N)
    x <- tryCatch(
      as.numeric(Matrix::solve(B, pik)),
      error = function(e) NULL
    )
    if (is.null(x) || length(x) != N || any(!is.finite(x))) next

    residual <- max(abs(as.numeric(B %*% x) - pik))
    if (min(x) < -feasibility_tol || residual > 10 * feasibility_tol) next

    x[abs(x) < feasibility_tol] <- 0
    x[x < 0 & x >= -feasibility_tol] <- 0

    residual2 <- max(abs(as.numeric(B %*% x) - pik))
    if (residual2 > 100 * feasibility_tol) next

    if (verbose) {
      message(sprintf(
        "Initial feasible basis found on attempt %d | perturbation %.3g | positive basic variables %d/%d",
        attempt, step, sum(x > feasibility_tol), N
      ))
    }

    return(list(
      samples = samples,
      x = x,
      B = B,
      attempt = attempt,
      perturbation = step,
      positive = sum(x > feasibility_tol),
      feasibility_error = residual2
    ))
  }

  stop(
    "Could not construct a nonsingular feasible initial basis. ",
    "Try increasing max_initial_tries or perturb_scale."
  )
}


# -----------------------------------------------------------------------------
# Product-form (eta) inverse for the revised simplex
# -----------------------------------------------------------------------------
# If column k of B is replaced by a and d = B^{-1}a, then
#
#       B_new = B E,
#
# where E is the identity matrix with column k replaced by d.  Between sparse
# LU refactorizations we therefore store only (k,d).  This avoids solving two
# newly factorized sparse systems at every simplex iteration.
#
# Periodic refactorization controls numerical drift.

.cos_factor_solve <- function(factor, A, rhs) {
  tryCatch(
    as.numeric(solve(factor, rhs)),
    error = function(e) as.numeric(Matrix::solve(A, rhs))
  )
}


.cos_linear_engine <- function(B) {
  state <- new.env(parent = emptyenv())

  state$refactor <- function(Bnew) {
    state$B0 <- methods::as(Bnew, "CsparseMatrix")
    state$B0T <- Matrix::t(state$B0)
    state$fac <- tryCatch(Matrix::lu(state$B0), error = function(e) NULL)
    state$facT <- tryCatch(Matrix::lu(state$B0T), error = function(e) NULL)
    state$etas <- list()
    invisible(NULL)
  }

  state$solve_base <- function(rhs) {
    if (!is.null(state$fac)) {
      .cos_factor_solve(state$fac, state$B0, rhs)
    } else {
      as.numeric(Matrix::solve(state$B0, rhs))
    }
  }

  state$solve_base_t <- function(rhs) {
    if (!is.null(state$facT)) {
      .cos_factor_solve(state$facT, state$B0T, rhs)
    } else {
      as.numeric(Matrix::solve(state$B0T, rhs))
    }
  }

  # Solve current B x = rhs.
  state$solve <- function(rhs) {
    z <- state$solve_base(rhs)
    if (length(state$etas)) {
      for (eta in state$etas) {
        k <- eta$k
        d <- eta$d
        if (abs(d[k]) < 1e-14) stop("Numerically singular eta update.")
        xk <- z[k] / d[k]
        z <- z - d * xk
        z[k] <- xk
      }
    }
    z
  }

  # Solve current B^T y = rhs.
  state$solve_t <- function(rhs) {
    z <- as.numeric(rhs)
    if (length(state$etas)) {
      for (ii in rev(seq_along(state$etas))) {
        eta <- state$etas[[ii]]
        k <- eta$k
        d <- eta$d
        if (abs(d[k]) < 1e-14) stop("Numerically singular transpose eta update.")

        other <- sum(d * z) - d[k] * z[k]
        z[k] <- (z[k] - other) / d[k]
      }
    }
    state$solve_base_t(z)
  }

  state$add_eta <- function(k, d) {
    state$etas[[length(state$etas) + 1L]] <- list(k = as.integer(k), d = as.numeric(d))
    invisible(NULL)
  }

  state$n_eta <- function() length(state$etas)

  state$refactor(B)
  state
}


# -----------------------------------------------------------------------------
# Pricing
# -----------------------------------------------------------------------------

.cos_price_heuristic <- function(
    z, dual, target, n, basis_keys,
    n_starts, max_swaps, improvement_tol) {

  ans <- cos_pricing_swap_cpp(
    z = z,
    dual = dual,
    target = target,
    n = as.integer(n),
    n_starts = as.integer(n_starts),
    max_swaps = as.integer(max_swaps),
    improvement_tol = improvement_tol
  )

  ord <- order(ans$reduced_costs)
  samples <- as.matrix(ans$samples)

  for (ii in ord) {
    s <- as.integer(samples[ii, ])
    key <- .cos_sample_key(s)
    if (!(key %in% basis_keys)) {
      return(list(
        sample = s,
        reduced_cost = as.numeric(ans$reduced_costs[ii]),
        certified = FALSE,
        method = "heuristic_swap",
        accepted_swaps = as.integer(ans$accepted_swaps[ii])
      ))
    }
  }

  list(
    sample = NULL,
    reduced_cost = Inf,
    certified = FALSE,
    method = "heuristic_swap",
    accepted_swaps = NA_integer_
  )
}


.cos_price_exact <- function(z, dual, target, n) {
  ans <- cos_pricing_exact_cpp(
    z = z,
    dual = dual,
    target = target,
    n = as.integer(n)
  )
  list(
    sample = as.integer(ans$sample),
    reduced_cost = as.numeric(ans$reduced_cost),
    certified = TRUE,
    method = "exact_enumeration",
    evaluated = ans$evaluated
  )
}


.cos_pricing <- function(
    pricing,
    z, dual, target, n, basis_keys,
    pricing_starts,
    pricing_max_swaps,
    pricing_tol,
    exact_threshold) {

  N <- length(z)
  exact_possible <- .cos_choose_within(N, n, exact_threshold)

  if (pricing == "exact") {
    if (!exact_possible) {
      stop(
        "Exact pricing would enumerate choose(N,n) samples, which exceeds exact_threshold. ",
        "Increase exact_threshold explicitly or use pricing='heuristic'/'hybrid'."
      )
    }
    return(.cos_price_exact(z, dual, target, n))
  }

  h <- .cos_price_heuristic(
    z, dual, target, n, basis_keys,
    n_starts = pricing_starts,
    max_swaps = pricing_max_swaps,
    improvement_tol = pricing_tol
  )

  # Any negative reduced-cost column is sufficient for a valid simplex pivot.
  if (!is.null(h$sample) && h$reduced_cost < -pricing_tol) return(h)

  if (pricing == "hybrid" && exact_possible) {
    return(.cos_price_exact(z, dual, target, n))
  }

  h
}


# -----------------------------------------------------------------------------
# Main revised-simplex optimizer
# -----------------------------------------------------------------------------

cos_design <- function(
    values,
    pik,
    pricing = c("hybrid", "heuristic", "exact"),
    max_iter = 500L,
    pricing_starts = 20L,
    pricing_max_swaps = 2000L,
    exact_threshold = 2e6,
    refactor_every = 25L,
    tol = 1e-9,
    pricing_tol = 1e-10,
    support_tol = 1e-10,
    pivot_tol = 1e-12,
    degenerate_limit = 100L,
    perturb_scale = 1e-7,
    max_initial_tries = 30L,
    seed = NULL,
    verbose = TRUE,
    label = "COS") {

  pricing <- match.arg(pricing)
  pr <- .cos_validate_problem(values, pik)
  values <- pr$values
  pik <- pr$pik
  N <- pr$N
  n <- pr$n

  max_iter <- as.integer(max_iter)
  refactor_every <- max(1L, as.integer(refactor_every))
  pricing_starts <- max(1L, as.integer(pricing_starts))
  pricing_max_swaps <- max(1L, as.integer(pricing_max_swaps))

  if (!is.null(seed)) set.seed(seed)

  # For the HT estimator, z_k = value_k / pi_k and the target is the total.
  z <- values / pik
  target <- sum(values)

  t0 <- proc.time()[["elapsed"]]

  init <- .cos_initial_basis(
    pik = pik,
    n = n,
    perturb_scale = perturb_scale,
    max_tries = max_initial_tries,
    feasibility_tol = max(tol, 1e-10),
    seed = seed,
    verbose = verbose
  )

  basis_samples <- init$samples
  x <- as.numeric(init$x)
  cB <- .cos_sample_costs(basis_samples, z, target)
  basis_keys <- apply(basis_samples, 1L, .cos_sample_key)

  engine <- .cos_linear_engine(init$B)

  trace_rows <- vector("list", max_iter + 1L)
  trace_n <- 0L
  status <- "max_iter"
  certified <- FALSE
  degenerate_count <- 0L

  if (verbose) {
    message(sprintf(
      "%s revised simplex | N=%d n=%d | pricing=%s | initial objective=%.8g",
      label, N, n, pricing, sum(x * cB)
    ))
  }

  for (iter in seq_len(max_iter)) {
    objective_before <- sum(x * cB)

    # Dual multipliers: B^T y = c_B.
    dual <- engine$solve_t(cB)

    price <- .cos_pricing(
      pricing = pricing,
      z = z,
      dual = dual,
      target = target,
      n = n,
      basis_keys = basis_keys,
      pricing_starts = pricing_starts,
      pricing_max_swaps = pricing_max_swaps,
      pricing_tol = pricing_tol,
      exact_threshold = exact_threshold
    )

    trace_n <- trace_n + 1L
    trace_rows[[trace_n]] <- data.frame(
      iteration = iter,
      objective = objective_before,
      reduced_cost = price$reduced_cost,
      pricing_method = price$method,
      positive_support = sum(x > support_tol),
      eta_updates = engine$n_eta(),
      elapsed_seconds = proc.time()[["elapsed"]] - t0,
      stringsAsFactors = FALSE
    )

    if (verbose && (iter == 1L || iter %% 10L == 0L || price$reduced_cost >= -pricing_tol)) {
      message(sprintf(
        "iter %d | objective %.8g | reduced cost %.4g | support %d",
        iter, objective_before, price$reduced_cost, sum(x > support_tol)
      ))
    }

    # Termination.
    if (is.null(price$sample) || !is.finite(price$reduced_cost) ||
        price$reduced_cost >= -pricing_tol) {
      if (isTRUE(price$certified)) {
        status <- "optimal_exact"
        certified <- TRUE
      } else {
        status <- "heuristic_stationary"
        certified <- FALSE
      }
      break
    }

    enter_sample <- sort(as.integer(price$sample))
    enter_key <- .cos_sample_key(enter_sample)

    # A genuinely negative reduced-cost basic column should be impossible except
    # for numerical error.  Refactor once before failing.
    if (enter_key %in% basis_keys) {
      Bnow <- .cos_basis_matrix(basis_samples, N)
      engine$refactor(Bnow)
      dual <- engine$solve_t(cB)
      rc_check <- (sum(z[enter_sample]) - target)^2 - sum(dual[enter_sample])
      if (rc_check >= -10 * pricing_tol) {
        status <- if (pricing == "exact") "optimal_exact" else "heuristic_stationary"
        certified <- identical(status, "optimal_exact")
        break
      }
      stop("Pricing returned a basic column with a materially negative reduced cost.")
    }

    a <- numeric(N)
    a[enter_sample] <- 1

    # Direction d = B^{-1} a.
    d <- engine$solve(a)
    eligible <- which(d > pivot_tol)
    if (!length(eligible)) {
      stop("No positive pivot direction found; numerical failure in the basis solve.")
    }

    ratios <- x[eligible] / d[eligible]
    min_ratio <- min(ratios)
    tie <- eligible[which(abs(ratios - min_ratio) <= max(tol, 1e-12))]
    leave <- min(tie)  # deterministic tie-breaking
    theta <- x[leave] / d[leave]

    if (!is.finite(theta) || theta < -100 * tol) {
      stop("Invalid ratio-test step; current basis appears numerically unstable.")
    }
    theta <- max(theta, 0)

    old_key <- basis_keys[leave]

    # New basic coordinates after replacing column 'leave' by the entering one.
    x <- x - theta * d
    x[abs(x) < support_tol] <- 0
    x[x < 0 & x > -100 * tol] <- 0
    x[leave] <- theta

    basis_samples[leave, ] <- enter_sample
    basis_keys[leave] <- enter_key
    cB[leave] <- (sum(z[enter_sample]) - target)^2

    if (theta <= support_tol) {
      degenerate_count <- degenerate_count + 1L
    } else {
      degenerate_count <- 0L
    }

    # Product-form inverse update, with periodic sparse refactorization.
    force_refactor <- (
      abs(d[leave]) < 1e-10 ||
      engine$n_eta() + 1L >= refactor_every
    )

    if (force_refactor) {
      Bnow <- .cos_basis_matrix(basis_samples, N)
      engine$refactor(Bnow)
    } else {
      engine$add_eta(leave, d)
    }

    if (degenerate_count >= degenerate_limit) {
      status <- "degenerate_stall"
      certified <- FALSE
      warning("Stopping after many consecutive degenerate pivots.")
      break
    }

    # Cheap internal consistency check from time to time.
    if (iter %% refactor_every == 0L) {
      Bcheck <- .cos_basis_matrix(basis_samples, N)
      feas_err <- max(abs(as.numeric(Bcheck %*% x) - pik))
      if (feas_err > 1e-6) {
        stop(sprintf("Feasibility drift became too large (%.3g).", feas_err))
      }
    }
  }

  if (status == "max_iter" && verbose) {
    message("Maximum number of simplex iterations reached.")
  }

  positive <- which(x > support_tol)
  support <- basis_samples[positive, , drop = FALSE]
  probabilities <- x[positive]
  objective <- sum(probabilities * .cos_sample_costs(support, z, target))

  # Final feasibility diagnostics are computed from the positive support only.
  inc <- numeric(N)
  for (j in seq_len(nrow(support))) {
    inc[support[j, ]] <- inc[support[j, ]] + probabilities[j]
  }

  trace <- if (trace_n) {
    do.call(rbind, trace_rows[seq_len(trace_n)])
  } else {
    data.frame()
  }

  ans <- list(
    call = match.call(),
    label = label,
    N = N,
    n = n,
    pik = pik,
    design_values = values,
    target = target,
    support = support,
    probabilities = probabilities,
    objective = objective,
    rmse = sqrt(max(objective, 0)),
    relative_rmse = if (abs(target) > 0) sqrt(max(objective, 0)) / abs(target) else NA_real_,
    status = status,
    certified = certified,
    pricing = pricing,
    iterations = trace_n,
    trace = trace,
    inclusion_probabilities = inc,
    max_inclusion_error = max(abs(inc - pik)),
    probability_sum = sum(probabilities),
    elapsed_seconds = proc.time()[["elapsed"]] - t0,
    initialization = list(
      attempt = init$attempt,
      perturbation = init$perturbation,
      positive_basic_variables = init$positive,
      feasibility_error = init$feasibility_error
    ),
    control = list(
      max_iter = max_iter,
      pricing_starts = pricing_starts,
      pricing_max_swaps = pricing_max_swaps,
      exact_threshold = exact_threshold,
      refactor_every = refactor_every,
      tol = tol,
      pricing_tol = pricing_tol,
      support_tol = support_tol
    )
  )
  class(ans) <- "cos_design"
  ans
}


acos_design <- function(predicted_values, pik, ..., label = "ACOS") {
  cos_design(
    values = predicted_values,
    pik = pik,
    ...,
    label = label
  )
}


# -----------------------------------------------------------------------------
# Evaluation and sampling helpers
# -----------------------------------------------------------------------------

evaluate_cos_design <- function(fit, values) {
  if (!inherits(fit, "cos_design")) stop("fit must be a cos_design object.")
  values <- as.numeric(values)
  if (length(values) != fit$N) stop("values must have length fit$N.")
  if (any(!is.finite(values))) stop("values must be finite.")

  z <- values / fit$pik
  target <- sum(values)
  costs <- .cos_sample_costs(fit$support, z, target)
  mse <- sum(fit$probabilities * costs)

  list(
    mse = mse,
    rmse = sqrt(max(mse, 0)),
    relative_rmse = if (abs(target) > 0) sqrt(max(mse, 0)) / abs(target) else NA_real_,
    target = target
  )
}


cos_diagnostics <- function(fit) {
  if (!inherits(fit, "cos_design")) stop("fit must be a cos_design object.")
  data.frame(
    N = fit$N,
    n = fit$n,
    support_size = nrow(fit$support),
    probability_sum = fit$probability_sum,
    max_inclusion_error = fit$max_inclusion_error,
    objective = fit$objective,
    rmse = fit$rmse,
    relative_rmse = fit$relative_rmse,
    iterations = fit$iterations,
    elapsed_seconds = fit$elapsed_seconds,
    status = fit$status,
    certified = fit$certified,
    pricing = fit$pricing,
    stringsAsFactors = FALSE
  )
}


draw_cos_sample <- function(fit, seed = NULL) {
  if (!inherits(fit, "cos_design")) stop("fit must be a cos_design object.")
  if (!is.null(seed)) set.seed(seed)
  p <- pmax(fit$probabilities, 0)
  p <- p / sum(p)
  j <- sample.int(nrow(fit$support), size = 1L, prob = p)
  as.integer(fit$support[j, ])
}


print.cos_design <- function(x, ...) {
  cat(sprintf("%s design\n", x$label))
  cat("----------------------------------------\n")
  cat(sprintf("N = %d | n = %d | support = %d\n", x$N, x$n, nrow(x$support)))
  cat(sprintf("objective (MSE): %.10g\n", x$objective))
  cat(sprintf("RMSE: %.10g\n", x$rmse))
  if (is.finite(x$relative_rmse)) {
    cat(sprintf("relative RMSE: %.6g\n", x$relative_rmse))
  }
  cat(sprintf("iterations: %d | elapsed: %.3f s\n", x$iterations, x$elapsed_seconds))
  cat(sprintf("status: %s | certified: %s\n", x$status, x$certified))
  cat(sprintf("max |pi_hat-pi|: %.3g | sum(p): %.12g\n",
              x$max_inclusion_error, x$probability_sum))
  invisible(x)
}


summary.cos_design <- function(object, ...) {
  cos_diagnostics(object)
}


# -----------------------------------------------------------------------------
# Small utilities for examples
# -----------------------------------------------------------------------------

cos_pps_inclusion_probabilities <- function(size, n, eps = 1e-10) {
  size <- as.numeric(size)
  N <- length(size)
  if (length(n) != 1L || n <= 0 || n >= N) stop("n must satisfy 0 < n < N.")
  if (any(!is.finite(size)) || any(size <= 0)) {
    stop("size must contain finite strictly positive values.")
  }

  pik <- numeric(N)
  active <- rep(TRUE, N)
  remaining <- as.numeric(n)

  repeat {
    idx <- which(active)
    if (!length(idx)) break

    cfac <- remaining / sum(size[idx])
    cand <- cfac * size[idx]
    certainty <- cand >= 1 - eps

    if (!any(certainty)) {
      pik[idx] <- cand
      break
    }

    cert_idx <- idx[certainty]
    pik[cert_idx] <- 1
    active[cert_idx] <- FALSE
    remaining <- remaining - length(cert_idx)

    if (remaining <= 0) break
  }

  # COS code assumes open probabilities; certainty units can simply be removed
  # from the optimization problem.  For the examples we reject that case.
  if (any(pik >= 1 - eps) || any(pik <= eps)) {
    stop("The requested PPS probabilities contain certainty/near-zero units.")
  }

  pik
}


cos_correlated_proxy <- function(y, rho, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (rho <= 0 || rho >= 1) stop("rho must lie in (0,1).")

  y <- as.numeric(y)
  yc <- y - mean(y)
  sy <- stats::sd(yc)
  noise <- rnorm(length(y))
  noise <- noise - mean(noise)

  # Orthogonalize the noise against y so the finite-population correlation is
  # controlled much more tightly than by repeated random trials.
  denom <- sum(yc * yc)
  if (denom <= 0) stop("y must have positive variance.")
  noise <- noise - sum(noise * yc) / denom * yc
  noise <- noise / stats::sd(noise)

  x_std <- rho * (yc / sy) + sqrt(1 - rho^2) * noise
  x <- mean(y) + sy * x_std
  x
}
