% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussianProcessSimulation.R
\name{simulateFunction}
\alias{simulateFunction}
\title{Simulation-based Function Generator}
\usage{
simulateFunction(
  object,
  nsim = 1,
  seed = NA,
  method = "spectral",
  xsim = NA,
  Ncos = 10,
  conditionalSimulation = TRUE
)
}
\arguments{
\item{object}{an object generated by \code{\link{gaussianProcessRegression}}}

\item{nsim}{the number of simulations, or test functions, to be created}

\item{seed}{a random number generator seed. Defaults to NA; which means no seed is set. For sake of reproducibility, set this to some integer value.}

\item{method}{\code{"decompose"} (default) or \code{"spectral"}, specifying the method used for simulation. 
Note that \code{"decompose"} is can be preferable, since it is exact but may be computationally infeasible for high-dimensional xsim.
On the other hand, \code{"spectral"} yields a function that can be evaluated at arbitrary sample locations.}

\item{xsim}{list of samples in input space, for simulation (only used for decomposition-based simulation, not for spectral method)}

\item{Ncos}{number of cosine functions (used with \code{method="spectral"} only)}

\item{conditionalSimulation}{whether (TRUE) or not (FALSE) to use conditional simulation}
}
\value{
a list of functions, where each function is the interpolation of one simulation realization. The length of the list depends on the nsim parameter.
}
\description{
Generate functions via simulation of GPR models, e.g.,
for assessment of optimization algorithms with
non-conditional or conditional simulation, based on real-world data.
}
\references{
N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.

C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
}
\seealso{
\code{\link{gaussianProcessRegression}}, \code{\link{simulate.cobbsGPR}}
}
