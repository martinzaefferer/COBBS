% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logging.R
\name{logfill}
\alias{logfill}
\title{Fill Logger}
\usage{
logfill(logs, freq)
}
\arguments{
\item{logs}{the data frame produced by repeated experiments with, using the logger wrapper}

\item{freq}{freq parameter of the logger}
}
\value{
a vector of scaled values
}
\description{
Fills logger result by last observation carried forward, to ensure that all runs (for each seed) have the same number of observations.
}
\examples{
## first experiment to be logged
lower <- -1
upper <- 3
fun <- function(x)sqrt(abs(x-1))
require(nloptr)
lbfgs(x0=-1,
      fn = logger(fun,"mytest",1),
      lower=lower,
      upper=upper)
## store result
log <- as.data.frame(get(paste("logger","y","mytest",sep=".")))
colnames(log) <- c("it","y","cy")
log$seed <- 1 
## Repeat with different initial guess
lower <- -1
upper <- 3
fun <- function(x)sqrt(abs(x-1))
lbfgs(x0=0,
      fn = logger(fun,"mytest",1),
      lower=lower,
      upper=upper)
## store result
log2 <- as.data.frame(get(paste("logger","y","mytest",sep=".")))
colnames(log2) <- c("it","y","cy")
log2$seed <- 2
## combine both results
log <- rbind(log,log2)
plot(log$it) # second run was longer
## fill
logf <- logfill(log,1)
plot(logf$it) # filled in values for first run

}
\keyword{internal}
