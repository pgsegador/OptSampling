# Multiple supplied starts, with optional continuation budgets per start.
vns_multistart <- function(pik, X, starts = NULL, n_starts = 10L,
                           time_budgets = Inf, seed = NULL, verbose = TRUE, ...) {
  pik <- .validate_pik(pik)
  N <- length(pik)
  X <- .validate_X(X, N)
  dots <- list(...)
  reserved <- c("start", "n_starts", "seed", "max_time", "verbose", "stratification")
  if (length(intersect(names(dots), reserved)))
    stop("Use starts, n_starts, time_budgets and seed in vns_multistart(); stratification is available in vns_systematic().")
  if (!is.numeric(time_budgets) || !length(time_budgets) ||
      anyNA(time_budgets) || any(time_budgets <= 0))
    stop("time_budgets must contain positive seconds (or Inf).")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(starts)) {
    if (length(n_starts) != 1L || !is.finite(n_starts) ||
        n_starts < 1 || n_starts != floor(n_starts))
      stop("n_starts must be a positive integer.")
    starts <- replicate(n_starts, .random_permutation(N), simplify = FALSE)
  } else {
    if (is.matrix(starts)) starts <- lapply(seq_len(nrow(starts)), function(i) starts[i, ])
    if (is.numeric(starts)) starts <- list(starts)
    if (!is.list(starts) || !length(starts)) stop("starts must contain at least one permutation.")
    starts <- lapply(starts, .validate_permutation, N = N)
  }
  fits <- vector("list", length(starts))
  attempts <- list()
  for (i in seq_along(starts)) {
    current <- starts[[i]]
    initial <- NULL
    for (b in seq_along(time_budgets)) {
      fit <- do.call(vns_systematic, c(list(pik = pik, X = X, start = current,
        max_time = time_budgets[b], verbose = FALSE), dots))
      if (is.null(initial)) initial <- fit$initial_objective
      d <- fit$diagnostics[fit$selected_runs, , drop = FALSE]
      attempts[[length(attempts) + 1L]] <- cbind(start_id = i, budget_id = b,
                                               time_budget = time_budgets[b], d)
      current <- fit$sigma
      if (isTRUE(d$type2_certified)) break
    }
    fit$initial_objective <- initial
    fit$improvement <- initial - fit$final_objective
    fit$relative_improvement <- if (abs(initial) > 1e-12) fit$improvement / abs(initial) else NA_real_
    fits[[i]] <- fit
    if (verbose) cat(sprintf("Start %d/%d: %.8g (%s)\n", i, length(starts), fit$final_objective, d$stop_reason))
  }
  best <- which.min(vapply(fits, function(x) x$final_objective, numeric(1)))
  out <- fits[[best]]
  out$call <- match.call()
  out$multistart <- list(best_start = best, fits = fits, starts = starts,
                         attempts = do.call(rbind, attempts), time_budgets = time_budgets)
  out
}

systematic_support <- function(pik, sigma = seq_along(pik), tol_pik = 1e-12) {
  pik <- .validate_pik(pik)
  sigma <- .validate_permutation(sigma, length(pik))
  systematic_support_cpp(sigma, pik, tol_pik)
}

draw_systematic_sample <- function(fit, seed = NULL) {
  if (!inherits(fit, "vns_systematic_fit")) stop("fit must be a vns_systematic_fit object.")
  if (!is.null(seed)) set.seed(seed)
  support <- systematic_support(fit$pik, fit$sigma)
  j <- sample.int(length(support$samples), 1L, prob = support$probabilities)
  as.integer(support$samples[[j]])
}
