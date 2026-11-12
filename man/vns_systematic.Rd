\name{vns_systematic}
\alias{vns_systematic}
\title{Optimize a Systematic Sampling Order}
\description{Optimize a Systematic Sampling Order.}
\usage{
vns_systematic (pik, X, objective = c("balance", "spread"), objective_control = list(), start = NULL, 
    search = c("adaptive", "type1", "type2"), stratification = FALSE, n_strata = NULL, 
    strata = NULL, strata_control = list(), n_starts = 1L, mixture_top_k = 1L, improvement_tol = 1e-10, 
    tol_pik = 1e-12, max_iter = 100000L, max_time = Inf, seed = NULL, verbose = TRUE) 
}
\arguments{
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{X}{Numeric matrix with one row per population unit.}
\item{objective}{Balance or spatial spreading loss.}
\item{objective_control}{Loss controls: metric and scale for balance; coordinates, metric and W for spreading.}
\item{start}{Initial permutation or systematic mixture.}
\item{search}{Adaptive, adjacent (type1), or all-transposition (type2) local search.}
\item{stratification}{Use locally optimized strata to construct an initial order.}
\item{n_strata}{Number of automatic warm-start strata; ignored when strata is supplied.}
\item{strata}{Vector of stratum labels, one per unit.}
\item{strata_control}{List with max_iter (default 500) and max_time (default 1 second) for each local warm-start search.}
\item{n_starts}{Number of random starts.}
\item{mixture_top_k}{Number of best mixture components to optimize as warm starts.}
\item{improvement_tol}{Absolute minimum improvement required for an accepted move.}
\item{tol_pik}{Numerical tolerance for inclusion probabilities and support breakpoints.}
\item{max_iter}{Maximum number of accepted moves or simplex iterations per call.}
\item{max_time}{Time limit in seconds per run.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A vns_systematic_fit object containing sigma, objective values, diagnostics and runs.}

\examples{
set.seed(4)
pik <- rep(0.25, 12)
X <- matrix(runif(24, 1, 4), ncol = 2)
fit <- vns_systematic(pik, X, seed = 5, verbose = FALSE)
summary(fit)
draw_systematic_sample(fit, seed = 6)
}
