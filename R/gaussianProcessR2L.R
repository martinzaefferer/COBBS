

###################################################################################
#' Build Gaussian Process 2-Level (GPR2L) Model
#'
#' This function builds ...
#' 
#' @param x design matrix (sample locations)
#' @param y vector of observations at \code{x}
#' @param control (list), with the options for the model building procedure.
#'
#' @return an object of class \code{cobbsGPRTE}.
#'
#' @export
#' @seealso \code{\link{predict.cobbsGPRTE}}, \code{\link{gaussianProcessRegression}}
###################################################################################
gaussianProcessR2L <- function(x,y,control=list()){
  npar <- ncol(x) #dimensionality
  con<-list(
    thetaLower=1e-6, 
    thetaUpper=1e12, 
    useLambda=FALSE,
    modelControl1=list(), #controls passed to the L1 model
    modelControl2=list(), #controls passed to the L2 model
    minsize=npar*20 #minimum number of observations for each partition to be modeled
  )
  con[names(control)] <- control
  control<-con
  
  ## TODO: defaults L1, useLambda, reinterpolate
  
  #control$modelControl$target <- c("y","s") #required for weighting the ensemble predictions
  
  control$modelControl1$useLambda <- TRUE
  control$modelControl1$thetaUpper <- control$thetaUpper
  control$modelControl1$thetaLower <- control$thetaLower
  control$modelControl1$lambdaUpper <- 0
  control$modelControl1$lambdaLower <- -1
  
  
  ## build global trend model 1
  model1 <- gaussianProcessRegression(x=x,y=y,control=control$modelControl1)
  
  
  control$modelControl2$useLambda <- control$useLambda
  control$modelControl2$thetaUpper <- control$thetaUpper
  control$modelControl2$thetaLower <- model1$dmodeltheta
  
  y2 <- y-predict(model1,x)$y
  
  ## build localized model 2
  model2 <- gaussianProcessRegression(x=x,y=y2,control=control$modelControl2)
  
  ## return
  fit <- list(fit1=model1,fit2=model2)
  class(fit)<- "cobbsGPR2L"
  return(fit)
}



###################################################################################
#' Predict cobbsGPRTE Model
#' 
#' Predict with model produced by \code{\link{gaussianProcessR2L}}.
#'
#' @param object GPRTE model (settings and parameters) of class \code{cobbsGPRTE}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted mean \code{y}, uncertainty / standard deviation \code{s} (optional) and expected improvement \code{ei} (optional). 
#' Whether \code{s} and \code{ei} are returned is specified by the vector of strings \code{object$target}, which then contains \code{"s"} and \code{"ei"}.
#'
#'
#' @seealso \code{\link{gaussianProcessRegression}}, \code{\link{gaussianProcessR2L}}
#' @export
#' @keywords internal
###################################################################################
predict.cobbsGPR2L <- function(object,newdata,...){
  
  ## predict with each model
  prediction1 <- predict(object$fit1,newdata)
  prediction2 <- predict(object$fit2,newdata)
  
  ## combine
  ret <- list(y=prediction1$y + prediction2$y)
  
  ## uncertainty
  if(!is.null(prediction1$s)&!is.null(prediction2$s))
    ret$s <- prediction1$s + prediction2$s

  return(ret)
}

###################################################################################
#' Gaussian Process 2-Level Simulation
#' 
#' Simulation of a 2-Level Gaussian process model.
#'
#' @param object GPR2L model (settings and parameters) of class \code{cobbsGPR2L}.
#' @param xsim list of samples in input space, to be simulated at
#' @param method \code{"decompose"} (default) or \code{"spectral"}, specifying the method used for simulation. 
#' Note that \code{"decompose"} is can be preferable, since it is exact but may be computationally infeasible for high-dimensional xsim.
#' On the other hand, \code{"spectral"} yields a function that can be evaluated at arbitrary sample locations.
#' @param nsim number of simulations
#' @param seed random number generator seed. Defaults to NA, in which case no seed is set
#' @param conditionalSimulation  logical, if set to TRUE (default), the simulation is conditioned with the training data of the GPR model.
#' Else, the simulation is non-conditional.
#' @param Ncos number of cosine functions (used with \code{method="spectral"} only)
#' @param returnAll if set to TRUE, a list with the simulated values (y) and the corresponding covariance matrix (covar)
#' of the simulated samples is returned. 
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$simulationReturnAll}
#'
#' @seealso \code{\link{gaussianProcessR2L}}, \code{\link{predict.cobbsGPR2L}}, \code{\link{simulate.cobbsGPR}}
#' @export
###################################################################################
simulate.cobbsGPR2L <- function(object,nsim=1,seed=NA,xsim=NA,method="spectral",conditionalSimulation=TRUE,Ncos=10,returnAll=FALSE,...){
  if(method=="decompose"){
    sim1 <- simulate(object$fit1,nsim=nsim,seed=seed,xsim=xsim,method="decompose",
                     conditionalSimulation=conditionalSimulation,Ncos=Ncos,returnAll=returnAll)
    sim2 <- simulate(object$fit2,nsim=nsim,seed=seed,xsim=xsim,method="decompose",
                     conditionalSimulation=conditionalSimulation,Ncos=Ncos,returnAll=returnAll)
    ## combine simulations from each sub-model
    simresult <- sim1+sim2 # TODO test
    #plot(xsim,simresult[,1],type="l")
    ## TODO: returnall? for function-like approach? see simulateFunction...
    return(simresult)
  }else if(method=="spectral"){
    simresult <- NULL
    simfun <- NULL
    
    for(i in 1:nsim){
      fun1 <- simulationSpectral(object=object$fit1,conditionalSimulation=conditionalSimulation,Ncos=Ncos)
      fun2 <- simulationSpectral(object=object$fit2,conditionalSimulation=conditionalSimulation,Ncos=Ncos)
      e <- new.env()
      e$fun1 <- fun1
      e$fun2 <- fun2
      sfun <- function(newdata){
        fun1(newdata) + fun2(newdata)
      }
      environment(sfun) <- e
      if(!is.na(xsim[1]))
        sres <- sfun(xsim)
      else 
        sres <- NA
      simresult <- cbind(simresult,sres)
      simfun <- c(simfun,sfun)
    }
    if(returnAll){
      return(list(y=simresult,simfun=simfun))
    }else{
      return(simresult)
    }
  }else{
    stop("The specified method in simulate.cobbsGPR2L does not exist. Use 'decompose' or 'spectral'")
  }  
}
