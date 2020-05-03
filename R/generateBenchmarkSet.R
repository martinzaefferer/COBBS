###################################################################################
#' Generate Benchmark Sets for Continuous Optimization
#' 
#' Based on provided data and parameters, this function trains a GPR model,
#' then generates the corresponding simulations, and
#' returns them as a set of functions to the user.
#' For the sake of comparison, the predictor (estimation) of the GPR model
#' is also returned
#' 
#' @param x matrix of sample locations
#' @param y vector of observations at sample locations x
#' @param control (list), specifying the options for model training and simulation:
#' \describe{
#' \item{\code{modelControl}}{a list of controls for the GPR model, passed to \code{control} argument of \code{\link{buildKriging}}.}
#' \item{\code{method}}{a string specifiying the simulation method: \code{"spectral"} (default) or \code{"decompose"}.}
#' \item{\code{xsim}}{Locations for simulation, for decompose method only. Defaults to \code{NA}.}
#' \item{\code{nsim}}{Number of simulation instances to generate. Defaults to \code{1}.}
#' \item{\code{Ncos}}{Parameter N, defining the number of cosine functions superimposed with \code{method="spectral"}. Default is \code{100}.}
#' \item{\code{conditionalSimulation}}{Boolean specifying whether conditional simulation is used or not. Default is \code{FALSE}.}
#' \item{\code{seed}}{Integer specifying a random number generator seed. Defaults to \code{NA} (no seed is set).}
#' }
#'
#' @return 
#' \describe{
#' \item{\code{fit}}{object/fit produced by \code{\link{buildKriging}}}
#' \item{\code{estimation}}{a function, GPR predictor} 
#' \item{\code{simulation}}{a list of one or more test functions to be used for benchmarking (based on GPR simulation)}
#' }
#' 
#' @examples
#' ## generate some data
#' seed <- 1234
#' groundtruth <- function(x){
#'   x=matrix(x,,1)
#'   apply(x,1,smoof::makeDeflectedCorrugatedSpringFunction(1))
#' }
#' lower = 0
#' upper = 10
#' set.seed(seed)
#' x <- matrix(c(runif(8,lower,upper)),,1)
#' y <- matrix(groundtruth(x),,1)
#' ## specify some model configuration
#' mc  <- list(useLambda=FALSE,budgetAlgTheta=200,thetaUpper=1e4)
#' ## and some configuration details for the simulation
#' cntrl <- list(modelControl=mc,
#'               nsim=1,
#'               seed=seed,
#'               method="spectral",
#'               #method="decompose",xsim=matrix(c(runif(200,-32,32)),,1),
#'               conditionalSimulation=TRUE
#' )
#' ## generate model and functions
#' cobbsResult <- generateCOBBS(x,y,cntrl)
#' cobbsResult$fit
#' ## plot the functions
#' require(ggplot2)
#' xseq <- seq(from=lower,by=0.001,to=upper)
#' p <- ggplot()
#' ## plot trained model (predictor, estimatation)
#' df <- data.frame(x=xseq,y=cobbsResult$estimation(xseq))
#' p <- p + geom_line(aes(x=x,y=y),df,size=1,color="darkgray")
#' ## plot  true function 
#' df <- data.frame(x=xseq,y=groundtruth(xseq))
#' p <- p + geom_line(aes(x=x,y=y),df,size=1,linetype=2)
#' ## plot training data
#' df <- data.frame(x=x,y=y)
#' p <- p + geom_point(aes(x=x,y=y),df,size=2,fill="white",shape = 21)
#' ## plot simulation
#' df <- data.frame(x=xseq,y=cobbsResult$simulation[[1]](xseq))
#' p <- p + geom_line(aes(x=x,y=y),df,size=1,color="red")
#' ## some stuff for good looks
#' p <- p + theme_classic()
#' p <- p + xlab("x") + ylab("f(x)") 
#' print(p)
#' ## prepare an expression that will be run during the experiments
#' ## here: L-BFGS-B with random initial guess
#' require(nloptr)
#' expr <- expression(
#'   res <- lbfgs(x0=matrix(runif(1,lower,upper),1,1),fn = fnlog,
#'                lower=lower,upper=upper)
#' )
#' ## run the experiments, with logging
#' ## with each objective function produced by COBBS
#' resgt <- loggedExperiment(expr, groundtruth, 1:10)
#' reses <- loggedExperiment(expr, cobbsResult$estimation, 1:10)
#' ressi <- loggedExperiment(expr, cobbsResult$simulation[[1]], 1:10)
#' ## plot results
#' print(plotBenchmarkPerformance(list(resgt,reses)))
#' print(plotBenchmarkPerformance(list(resgt,reses,ressi),c("groundtruth","estimation","simulation")))
#' ## plot error, comparing against groundtruth
#' print(plotBenchmarkValidation(resgt,list(reses),c("estimation")))
#' print(plotBenchmarkValidation(resgt,list(reses,ressi),c("estimation","simulation")))
#'
#' @export
###################################################################################
generateCOBBS <- function(x,y,control=list()){
  ## default settings
  con<-list(
    modelControl = list(),
    method="spectral",
    nsim=1,
    xsim=NA,
    Ncos=100,
    conditionalSimulation=FALSE,
    seed=NA
  )
  con[names(control)] <- control
  control<-con
  rm(con)
  
  fit <- buildKriging(x,y,control=control$modelControl)
  estimFun <- evaluateModel(fit)
  simFun <- simulateFunction(object=fit,
                             method=control$method, 
                             xsim=control$xsim, 
                             nsim=control$nsim, 
                             Ncos=control$Ncos,
                             conditionalSimulation=control$conditionalSimulation,
                             seed=control$seed)
  
  res <- list(fit=fit,estimation=estimFun,simulation=simFun)
  res #return
}
# 
# ## example plots 1dim
# lower = -1
# upper = 1
# expr <- expression(
#   res <- lbfgs(x0=runif(1,lower,upper),fn = fnlog,
#                      lower=lower,upper=upper)
# )
# res1 <- loggedExperiment(expr, function(x)sqrt(abs(x)), 1:10)
# res2 <- loggedExperiment(expr, function(x)x^2, 1:10)
# res3 <- loggedExperiment(expr, function(x)log(abs(x)+0.1), 1:10)
# 
# p <- plotBenchmarkValidation(res1,list(res2,res3),names=c("square","log"))
# print(p)
# p <- plotBenchmarkPerformance(list(res1,res2,res3),names=c("root","square","log"))
# print(p)
# 
# 
# ## example plots 2dim
# lower = c(-1,-1)
# upper = c(1,1)
# expr <- expression(
#   res <- lbfgs(x0=runif(2,lower,upper),fn = fnlog,
#                lower=lower,upper=upper)
# )
# res1 <- loggedExperiment(expr, function(x)sqrt(sum(abs(x))), 1:10)
# res2 <- loggedExperiment(expr, function(x)sum(x^2), 1:10)
# res3 <- loggedExperiment(expr, function(x)log(sum(abs(x)+0.1)), 1:10)
# 
# p <- plotBenchmarkValidation(res1,list(res2,res3),names=c("square","log"))
# print(p)
# p <- plotBenchmarkPerformance(list(res1,res2,res3),names=c("root","square","log"))
# print(p)

