

###################################################################################
#' Build Gaussian Process Regression Tree Ensemble (GPRTE) Model
#'
#' This function builds an ensemble of Gaussian Process model, where
#' each individual model is fitted to a partition of the parameter space.
#' Partitions are generated by a Tree-based approach, sequentially
#' dividing up the parameter space by making splits in single dimensions.
#' 
#' @param x design matrix (sample locations)
#' @param y vector of observations at \code{x}#' @param control (list), with the options for the model building procedure.
#'
#' @return an object of class \code{cobbsGPRTE}.
#'
#' @export
#' @seealso \code{\link{predict.cobbsGPRTE}}, \code{\link{gaussianProcessRegression}}
###################################################################################
gaussianProcessRTE <- function(x,y,control=list()){
  npar <- ncol(x) #dimensionality
  con<-list(
    modelControl=list(), #controls passed to the core model
    minsize=npar*20 #minimum number of observations for each partition to be modeled
  )
  con[names(control)] <- control
  control<-con
  
  control$modelControl$target <- c("y","s") #required for weighting the ensemble predictions
  
  ## split training data
  
  #####
  #d1 <- as.matrix(dist(x))
  #d2 <- as.matrix(dist(y))
  #dd <- d2/d1
  #vd <- apply(dd,1,var,na.rm=T)
  vd <- y
  df <- data.frame(x,vd)
  
  require(rpart)
  ## TODO: control the maximum node size, for complexity reduction purposes?
  fitTree <- rpart(vd~x,data=df,control=list(
    minsplit = 1,
    minbucket = control$minsize,
    maxdepth=30,
    cp=0
  ))
  ## we are not interested in regressing vd, rather predicting the leafe node
  ## hence:
  fitTree$frame$yval = order(as.numeric(rownames(fitTree$frame)))
  preds = predict(fitTree, newdata=df)  
  upred <- unique(preds)
  
  ## build a model in each partition
  models <- list()
  for(i in 1:length(upred)){
    selected <- preds == upred[i]
    xi <- x[selected,]
    yi <- y[selected,,drop=F]
    models[[i]] <-  gaussianProcessRegression(x=xi,y=yi,control=control$modelControl)
  }
  
  ## ...
  fit <- list(fits=models,fitTree=fitTree)
  class(fit)<- "cobbsGPRTE"
  return(fit)
}



###################################################################################
#' Predict cobbsGPRTE Model
#' 
#' Predict with potentially non-stationry ensemble model produced by \code{\link{gaussianProcessRTE}}.
#'
#' @param object GPRTE model (settings and parameters) of class \code{cobbsGPRTE}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted mean \code{y}, uncertainty / standard deviation \code{s} (optional) and expected improvement \code{ei} (optional). 
#' Whether \code{s} and \code{ei} are returned is specified by the vector of strings \code{object$target}, which then contains \code{"s"} and \code{"ei"}.
#'
#'
#' @seealso \code{\link{gaussianProcessRegression}}, \code{\link{gaussianProcessRTE}}
#' @export
#' @keywords internal
###################################################################################
predict.cobbsGPRTE <- function(object,newdata,...){
  
  ## predict with each model
  predictions <- list()
  for(i in 1:length(object$fits)){
    fit <- object$fits[[i]]
    predictions[[i]] <- predict(fit,newdata)
  }
  
  ps <- do.call(rbind,sapply(predictions,'[',"s"))
  psnegsquare <- ps^-2
  weights <- t(t(psnegsquare)/colSums(psnegsquare))
  py <- do.call(rbind,sapply(predictions,'[',"y"))
  #browser()
  ensembley <- colSums(py*weights)
  list(y=ensembley)
  #plot(ensembley)
  ## determine weights
}

###################################################################################
#' Gaussian Process Ensemble Simulation
#' 
#' Simulation of an ensemble of Gaussian process models.
#'
#' @param object GPRTE model (settings and parameters) of class \code{cobbsGPRTE}.
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
#' @seealso \code{\link{gaussianProcessRTE}}, \code{\link{predict.cobbsGPRTE}}, \code{\link{simulate.cobbsGPR}}
#' @export
###################################################################################
simulate.cobbsGPRTE <- function(object,nsim=1,seed=NA,xsim=NA,method="decompose",conditionalSimulation=TRUE,Ncos=10,returnAll=FALSE,...){
  if(method=="decompose"){
    ## predict with each sub-model, for weighting.
    predictions <- list()
    for(i in 1:length(object$fits)){
      fit <- object$fits[[i]]
      predictions[[i]] <- predict(fit,xsim) #TODO: xsim may be unavailable?
    }
    ## determine weights from predictions
    ps <- do.call(rbind,sapply(predictions,'[',"s"))
    psnegsquare <- ps^-2
    weights <- (t(psnegsquare)/colSums(psnegsquare))
    
    ## simulate each sub-model
    simulations <- list()
    for(i in 1:length(object$fits)){
      fit <- object$fits[[i]]
      simulations[[i]] <- simulate(object=fit,nsim=nsim,seed=seed,xsim=xsim,method=method,
                                   conditionalSimulation=conditionalSimulation,Ncos=Ncos,returnAll=returnAll)
      #plot(xsim,simulations[[i]][,1],type="l")
      simulations[[i]] <- simulations[[i]] * weights[,i] # include weights
    }  
    ## combine simulations from each sub-model
    simresult <- Reduce('+',simulations)
    #plot(xsim,simresult[,1],type="l")
    ## TODO: returnall? for function-like approach? see simulateFunction...
    return(simresult)
  }else if(method=="spectral"){
    simresult <- NULL
    simfun <- NULL
    for(i in 1:nsim){
      funs <- NULL
      for(i in 1:length(object$fits)){
        fit <- object$fits[[i]]
        res <- COBBS:::simulationSpectral(object=fit,conditionalSimulation=conditionalSimulation,Ncos=Ncos)
        funs <- c(funs,res)
      }
      fit
      funs
      object
      e <- new.env()
      e$funs <- funs
      e$object <- object
      sfun <- function(newdata){
        predictions <- list()
        for(i in 1:length(object$fits)){
          fit <- object$fits[[i]]
          predictions[[i]] <- predict(fit,newdata)
        }
        
        ps <- do.call(rbind,sapply(predictions,'[',"s"))
        psnegsquare <- ps^-2
        weights <- t(t(psnegsquare)/colSums(psnegsquare))
        py <- t(sapply(funs,function(x)x(newdata)))
        #browser()
        ensembley <- colSums(py*weights)
      }
      environment(sfun) <- e
      #
      #
      #
      #
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
    stop("The specified method in simulate.cobbsGPRTE does not exist. Use 'decompose' or 'spectral'")
  }  
}

## TODO: function&spectral simulation
## TODO: ensemble uncertainty
## TODO: truly noise data will confuse the tree building process. the "nonlinearity" measure is not robust to noise.
## TODO: if one of the models is extremely good, it may underestimate 
## the predicted variance even in unobserved regions. 
## then, the ensemble fails spectacularly
## e.g., this happens when some region has linear data, or worse, constant
## solution might be a fully bayesian treatment of the parameters, marginilization... with mcmc?
