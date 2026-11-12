\name{estimate_sampling_design_mc}
\alias{estimate_sampling_design_mc}
\title{Estimate Sampling Design Moments by Simulation}
\description{Estimate Sampling Design Moments by Simulation.}
\usage{
estimate_sampling_design_mc (sampler, N, M = 2000L, n_pairs = 2000L, pairs = NULL, seed = NULL, verbose = TRUE) 
}
\arguments{
\item{sampler}{Function returning unit indices or a binary incidence vector.}
\item{N}{Population size.}
\item{M}{Number of Monte Carlo replications.}
\item{n_pairs}{Number of distinct population pairs used for moment estimation.}
\item{pairs}{Optional two-column matrix of unit-label pairs.}
\item{seed}{Optional random seed.}
\item{verbose}{Print progress when TRUE.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{A numeric result, design object, or diagnostics as described in the README.}

