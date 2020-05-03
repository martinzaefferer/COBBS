

###################################################################################
#' Function Evaluation Logger
#' 
#' Wraps your objective function, to Log each function evaluation during an optimization experiment.
#' 
#' @param fn the function to be monitored by the logger (an objective function, following the schema y=f(x)).
#' @param id an ID that will be part of the names of the variables that store the logging information
#' @param freq a positive integer, giving the frequency of the logging: the history of y (and possibly x) values is permanently recorded every freq evaluations. Note: the best observed value is updated in every evaluation, regardless of this value. Default: 1 (every iteration is stored by the logger).
#' @param logx boolean, specifying whether x values will be logged as well (default: FALSE).
#' @param envir the environment in which the logger variables are stored (default: .GlobalEnv)
#'
#' @return A function which takes the same arguments as fn, and returns the same values. In addition, it writes some variables to the specified environment.
#' 
#' @examples
#' lower <- -1
#' upper <- 3
#' fun <- function(x)sqrt(abs(x-1))
#' require(nloptr)
#' lbfgs(x0=-1,
#'       fn = logger(fun,"mytest",1),
#'       lower=lower,
#'       upper=upper)
#' ## see the logged minimum value:
#' logger.min.mytest
#' ## plot the progress over time
#' plot(logger.y.mytest[,c(1,2)])
#' plot(logger.y.mytest[,c(1,3)])
#'
#' @export
###################################################################################
logger <- function(fn,id,freq=1,logx=FALSE,envir=.GlobalEnv){
  force(fn)
  force(id)
  force(logx)
  if(logx){
    assign(paste("logger","y",id,sep="."),c(),envir=envir)
    assign(paste("logger","x",id,sep="."),NULL,envir=envir)
    assign(paste("logger","i",id,sep="."),0,envir=envir)      
    assign(paste("logger","min",id,sep="."),Inf,envir=envir) 
    resfun <- function(x){
      i <- get(paste("logger","i",id,sep="."),envir=envir)
      currentMinimum <- get(paste("logger","min",id,sep="."),envir=envir)
      y <- fn(x)
      if(y <= currentMinimum){
        currentMinimum <- y
        assign(paste("logger","min",id,sep="."),y,envir=envir) 
      }
      assign(paste("logger","i",id,sep="."),i+length(y),envir=envir)
      if(((i+1) %% freq)==0){
        yhist <- get(paste("logger","y",id,sep="."),envir=envir)
        xhist <- get(paste("logger","x",id,sep="."),envir=envir)
        xhist <- rbind(xhist,x)
        yhist <- rbind(yhist,c(i+1,y,currentMinimum))
        assign(paste("logger","y",id,sep="."),yhist,envir=envir)
        assign(paste("logger","x",id,sep="."),xhist,envir=envir)
      }
      y
    }    
  }else{
    assign(paste("logger","y",id,sep="."),NULL,envir=envir)
    assign(paste("logger","i",id,sep="."),0,envir=envir) 
    assign(paste("logger","min",id,sep="."),Inf,envir=envir) 
    resfun <- function(x){
      i <- get(paste("logger","i",id,sep="."),envir=envir)
      currentMinimum <- get(paste("logger","min",id,sep="."),envir=envir)
      y <- fn(x)
      if(y <= currentMinimum){
        currentMinimum <- y
        assign(paste("logger","min",id,sep="."),y,envir=envir) 
      }
      assign(paste("logger","i",id,sep="."),i+length(y),envir=envir)
      if(((i+1) %% freq)==0){
        yhist <- get(paste("logger","y",id,sep="."),envir=envir)
        yhist <- rbind(yhist,c(i+1,y,currentMinimum))
        assign(paste("logger","y",id,sep="."),yhist,envir=envir)
      }
      y
    }
  }
  force(resfun)
  resfun
}

###################################################################################
#' Scaling to [0,1]
#' 
#' Scale a set of values to the interval from zero to one.
#' 
#' @param x a vector of values to be scaled
#'
#' @return a vector of scaled values
#' 
#' @examples
#' scale01(1:10)
#'
#' @export
#' @keywords internal
###################################################################################
scale01 <- function(x){
  x <- x-min(x)
  x <- x/max(x)
  x
}




###################################################################################
#' Fill Logger
#' 
#' Fills logger result by last observation carried forward, to ensure that all runs (for each seed) have the same number of observations.
#' 
#' @param logs the data frame produced by repeated experiments with, using the logger wrapper
#' @param freq freq parameter of the logger
#'
#' @return a vector of scaled values
#' 
#' @examples 
#' ## first experiment to be logged
#' lower <- -1
#' upper <- 3
#' fun <- function(x)sqrt(abs(x-1))
#' require(nloptr)
#' lbfgs(x0=-1,
#'       fn = logger(fun,"mytest",1),
#'       lower=lower,
#'       upper=upper)
#' ## store result
#' log <- as.data.frame(get(paste("logger","y","mytest",sep=".")))
#' colnames(log) <- c("it","y","cy")
#' log$seed <- 1 
#' ## Repeat with different initial guess
#' lower <- -1
#' upper <- 3
#' fun <- function(x)sqrt(abs(x-1))
#' lbfgs(x0=0,
#'       fn = logger(fun,"mytest",1),
#'       lower=lower,
#'       upper=upper)
#' ## store result
#' log2 <- as.data.frame(get(paste("logger","y","mytest",sep=".")))
#' colnames(log2) <- c("it","y","cy")
#' log2$seed <- 2
#' ## combine both results
#' log <- rbind(log,log2)
#' plot(log$it) # second run was longer
#' ## fill
#' logf <- logfill(log,1)
#' plot(logf$it) # filled in values for first run
#'
#' @export
#' @keywords internal
###################################################################################
logfill <- function(logs,freq){
  maxit <- max(logs$it)
  logres <- NULL
  for(seed in unique(logs$seed)){
    log <- logs[logs$seed == seed,]
    if(tail(log$it,1)<maxit){
      its <- seq(from=freq,by=freq,to=maxit)
      its <- its[!(its %in% log$it)]
      logtmp <- cbind(data.frame(it=its,y=NA,cy=min(log$cy),seed=seed),matrix(NA,1,ncol(logs)-4))
      names(logtmp) <- names(log)
      log <- rbind(log,logtmp)
    }
    logres <- rbind(logres,log)
  }
  logres
}


###################################################################################
#' Logged Experimental Runs
#' 
#' This function performs optimization experiments (described by the expression expr),
#' where the objective function evaluations are logged as prescribed by the given arguments.
#' 
#' @param expr this R expression will be run (with \code{eval}), to perform a single optimization experiment.
#' Importantly, the expression should contain a call to a solver, which includes as an argument the objective function,
#' with the exact name \code{fnlog}. This function will be wrapped by \code{\link{logger}}. See the examples.
#' @param fun the actual objective function (\code{class(fun) = function}). 
#' @param seeds the random number seeds to be set before each experiment. expr will be called as often as length(seeds). Each call is preceded by a call to set.seed, with the respective element of this vector.
#' @param freq a positive integer, giving the frequency of the logging: the history of y (and possibly x) values is permanently recorded every freq evaluations. Note: the best observed value is updated in every evaluation, regardless of this value. Default: 1 (every iteration is stored by the logger).
#' @param logx boolean, specifying whether x values will be logged as well (default: FALSE).
#'
#' @return a data.frame containing the logged results: column 'it' gives iteration numbers, 'y' gives observations, 'cy' is the cumulated minimum of the observations over iterations, and 'seed' gives the rng seed of the corresponding experiment.
#' 
#' @examples
#' lower <- -1
#' upper <- 1
#' require(nloptr)
#' expr <- expression(
#'   res <- lbfgs(x0=runif(1,lower,upper),
#'                fn = fnlog, #important: use this same variable name 'fnlog'
#'                lower=lower,upper=upper)
#' )
#' res <- loggedExperiment(expr, function(x)sqrt(abs(x)), 1:20)
#' plot(res[,c(1,3)],col=res$seed)
#'
#' @export
###################################################################################
loggedExperiment <- function(expr, fun, seeds, freq=1,logx=FALSE){
  logall <- NULL
  for(seed in seeds){
    set.seed(seed)
    fnlog <- logger(fun,paste("test",seed,sep="."),freq,logx,envir = environment())
    eval(expr,envir = environment())
    log <- as.data.frame(get(paste("logger","y","test",seed,sep="."),envir=environment()))
    if(logx){
      logxvals <- as.data.frame(get(paste("logger","x","test",seed,sep="."),envir=environment()))
      log <- cbind(log,logxvals)
      colnames(log) <- c("it","y","cy",paste("x",1:ncol(logxvals),sep=""))
    }else{
      colnames(log) <- c("it","y","cy")
    }
    log$seed <- seed
    logall <- rbind(logall,log)
  }
  logall <- logfill(logall,freq)
  return(logall)
}
