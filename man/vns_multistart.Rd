\name{vns_multistart}
\alias{vns_multistart}
\title{Optimize Several Starts with Continuation Budgets}
\description{Optimize Several Starts with Continuation Budgets.}
\usage{
vns_multistart (pik, X, starts = NULL, n_starts = 10L, time_budgets = Inf, seed = NULL, verbose = TRUE, 
    ...) 
}
\arguments{
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{X}{Numeric matrix with one row per population unit.}
\item{starts}{List of permutations, matrix with one permutation per row, or NULL for random starts.}
\item{n_starts}{Number of random starts.}
\item{time_budgets}{Positive time limits in seconds per continuation, not cumulative; Inf is allowed.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
\item{...}{Additional arguments passed to the underlying optimizer.}
}
\details{Continues each supplied or generated start using each time budget in sequence, stopping that start on a Type-II certificate. Chooses the smallest final objective. The multistart field retains starts, final fits, and all attempt diagnostics. Execution is sequential. Explicit starts determine the number of runs.}
\value{A vns_systematic_fit object containing sigma, objective values, diagnostics and runs.}

