\name{cos_design}
\alias{cos_design}
\title{Construct a Conditional Optimal Sampling Design}
\description{Construct a Conditional Optimal Sampling Design.}
\usage{
cos_design (values, pik, pricing = c("hybrid", "heuristic", "exact"), max_iter = 500L, pricing_starts = 20L, 
    pricing_max_swaps = 2000L, exact_threshold = 2e+06, refactor_every = 25L, tol = 1e-09, 
    pricing_tol = 1e-10, support_tol = 1e-10, pivot_tol = 1e-12, degenerate_limit = 100L, 
    perturb_scale = 1e-07, max_initial_tries = 30L, seed = NULL, verbose = TRUE, label = "COS") 
}
\arguments{
\item{values}{Finite population values for the total to be estimated.}
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{pricing}{Exact enumeration, heuristic swap pricing, or hybrid selection.}
\item{max_iter}{Maximum number of accepted moves or simplex iterations per call.}
\item{pricing_starts}{Number of independent starts in heuristic column pricing.}
\item{pricing_max_swaps}{Maximum swaps per heuristic pricing start.}
\item{exact_threshold}{Maximum number of possible samples for exact pricing in hybrid mode.}
\item{refactor_every}{Maximum eta updates between full basis refactorizations.}
\item{tol}{Numerical feasibility and ratio-test tolerance.}
\item{pricing_tol}{Reduced-cost tolerance for accepting an entering sample.}
\item{support_tol}{Threshold for dropping negligible support probabilities.}
\item{pivot_tol}{Minimum positive pivot direction for the ratio test.}
\item{degenerate_limit}{Stop after this many consecutive degenerate pivots.}
\item{perturb_scale}{Initial perturbation magnitude for constructing a feasible basis.}
\item{max_initial_tries}{Maximum attempts to construct an initial feasible basis.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
\item{label}{Display label stored in the fitted design.}
}
\details{COS requires strictly interior inclusion probabilities. Exact pricing may enumerate exponentially many samples. Heuristic termination is not a global optimality certificate. The returned certified and status fields describe termination. Pricing controls set restart and swap counts; exact_threshold limits hybrid enumeration. Numerical controls govern feasibility, pivots, support truncation, refactorization and degeneracy.}
\value{A cos_design object containing support, probabilities, objective, status, certified and trace.}

