\name{tsp_order}
\alias{tsp_order}
\title{Construct a Traveling-Salesman Warm Start}
\description{Construct a Traveling-Salesman Warm Start.}
\usage{
tsp_order (D, n_starts = 20L, method = "two_opt", seed = NULL, verbose = TRUE) 
}
\arguments{
\item{D}{Distance matrix for the traveling-salesman warm start.}
\item{n_starts}{Number of random starts.}
\item{method}{Method passed to TSP::solve_TSP; defaults to two_opt.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A numeric result, design object, or diagnostics as described in the README.}

