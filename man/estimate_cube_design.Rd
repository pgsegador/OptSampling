\name{estimate_cube_design}
\alias{estimate_cube_design}
\title{Estimate Cube Design Moments by Simulation}
\description{Estimate Cube Design Moments by Simulation.}
\usage{
estimate_cube_design (pik, X_balance, M = 2000L, n_pairs = 2000L, pairs = NULL, force_fixed_size = TRUE, 
    seed = NULL, verbose = TRUE) 
}
\arguments{
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{X_balance}{Matrix of cube balancing variables, one row per population unit.}
\item{M}{Number of Monte Carlo replications.}
\item{n_pairs}{Number of distinct population pairs used for moment estimation.}
\item{pairs}{Optional two-column matrix of unit-label pairs.}
\item{force_fixed_size}{Add pik as a cube balancing variable to enforce fixed sample size.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A numeric result, design object, or diagnostics as described in the README.}

