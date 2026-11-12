\name{vns_objective}
\alias{vns_objective}
\title{Evaluate the Exact Expected Systematic Loss}
\description{Evaluate the Exact Expected Systematic Loss.}
\usage{
vns_objective (sigma, pik, X, objective = c("balance", "spread"), objective_control = list(), 
    tol_pik = 1e-12) 
}
\arguments{
\item{sigma}{Permutation of the population labels, from 1 to N.}
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{X}{Numeric matrix with one row per population unit.}
\item{objective}{Balance or spatial spreading loss.}
\item{objective_control}{Loss controls: metric and scale for balance; coordinates, metric and W for spreading.}
\item{tol_pik}{Numerical tolerance for inclusion probabilities and support breakpoints.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A numeric result, design object, or diagnostics as described in the README.}

