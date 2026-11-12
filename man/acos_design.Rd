\name{acos_design}
\alias{acos_design}
\title{Construct a Design from Predicted Values}
\description{Construct a Design from Predicted Values.}
\usage{
acos_design (predicted_values, pik, ..., label = "ACOS") 
}
\arguments{
\item{predicted_values}{Design-stage predictions of the study variable.}
\item{pik}{Prescribed first-order inclusion probabilities, summing to an integer.}
\item{...}{Additional arguments passed to the underlying optimizer.}
\item{label}{Display label stored in the fitted design.}
}
\details{COS requires strictly interior inclusion probabilities. Exact pricing may enumerate exponentially many samples. Heuristic termination is not a global optimality certificate. The returned certified and status fields describe termination. Pricing controls set restart and swap counts; exact_threshold limits hybrid enumeration. Numerical controls govern feasibility, pivots, support truncation, refactorization and degeneracy.}
\value{A cos_design object containing support, probabilities, objective, status, certified and trace.}

