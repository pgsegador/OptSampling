\name{systematic_support}
\alias{systematic_support}
\title{Enumerate a Systematic Design Support}
\description{Enumerate a Systematic Design Support.}
\usage{
systematic_support (pik, sigma = seq_along(pik), tol_pik = 1e-12) 
}
\arguments{
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{sigma}{Permutation of the population labels, from 1 to N.}
\item{tol_pik}{Numerical tolerance for inclusion probabilities and support breakpoints.}
}
\details{See the package README for the tutorial, metric definitions, output fields, optional dependencies and search controls. VNS certificates concern local minima over the selected neighborhood and numerical tolerance, not global optimality.}
\value{List with samples (a list of integer vectors) and probabilities.}

