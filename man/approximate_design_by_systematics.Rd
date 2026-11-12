\name{approximate_design_by_systematics}
\alias{approximate_design_by_systematics}
\title{Approximate Second-Order Moments by a Systematic Mixture}
\description{Approximate Second-Order Moments by a Systematic Mixture.}
\usage{
approximate_design_by_systematics (design, pik, n_candidates = 100L, candidate_permutations = NULL, max_components = 20L, 
    alpha_tol = 1e-08, regression_max_iter = 5000L, regression_tol = 1e-09, tol_pik = 1e-12, 
    seed = NULL, verbose = TRUE) 
}
\arguments{
\item{design}{Empirical sampling design object.}
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{n_candidates}{Target number of candidate systematic orderings.}
\item{candidate_permutations}{Optional list or matrix of supplied orders; each matrix row is a permutation.}
\item{max_components}{Maximum number of mixture components retained after fitting.}
\item{alpha_tol}{Threshold for retaining fitted mixture weights.}
\item{regression_max_iter}{Maximum projected-gradient iterations for the simplex-constrained least-squares fit.}
\item{regression_tol}{Convergence tolerance on successive mixture weights.}
\item{tol_pik}{Numerical tolerance for inclusion probabilities and support breakpoints.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A numeric result, design object, or diagnostics as described in the README.}

