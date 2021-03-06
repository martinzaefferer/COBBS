
###################################################################################
#' Build Gaussian Process Regression (GPR) Model
#'
#' This function builds a GPR model (aka Kriging) based on code by Forrester et al..
#' By default exponents (p) are fixed at a value of two, and a nugget (or regularization constant) is used.
#' To correct the uncertainty estimates in case of nugget, re-interpolation is also by default turned on.
#'
#' The model uses a Gaussian kernel: \code{k(x,z)=exp(-sum(theta_i * |x_i-z_i|^p_i))}. By default, \code{p_i = 2}.
#' Note that if dimension \code{x_i} is a factor variable (see parameter \code{types}), Hamming distance will be used 
#' instead of \code{|x_i-z_i|}.
#' 
#' @param x design matrix (sample locations)
#' @param y vector of observations at \code{x}
#' @param control (list), with the options for the model building procedure:\cr
#' \code{types} a character vector giving the data type of each variable. All but "factor" will be handled as numeric, "factor" (categorical) variables will be subject to the hamming distance.\cr
#' \code{thetaLower} lower boundary for theta, default is \code{1e-6}\cr
#' \code{thetaUpper} upper boundary for theta, default is \code{1e12}\cr
#' \code{algTheta}  algorithm used to find theta, default is \code{DEinterface}.\cr
#' \code{budgetAlgTheta} budget for the above mentioned algorithm, default is \code{200}. The value will be multiplied with the length of the model parameter vector to be optimized.\cr
#' \code{optimizeP} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. Default is \code{FALSE}\cr
#' \code{useLambda} whether or not to use the regularization constant lambda (nugget effect). Default is \code{TRUE}.\cr
#' \code{lambdaLower} lower boundary for log10{lambda}, default is \code{-6}\cr 
#' \code{lambdaUpper} upper boundary for log10{lambda}, default is \code{0}\cr
#' \code{startTheta} optional start value for theta optimization, default is \code{NULL}\cr
#' \code{reinterpolate} whether (\code{TRUE}) or not (\code{FALSE},default) reinterpolation should be performed.
#' \code{target} target values of the prediction, a vector of strings. Each string specifies a value to be predicted, e.g., "y" for mean, "s" for standard deviation, "ei" for expected improvement. See also \code{\link{predict.cobbsGPR}}.
#' This can also be changed after the model has been built, by manipulating the respective \code{object$target} value.
#'
#' @return an object of class \code{cobbsGPR}. Basically a list, with the options and found parameters for the model which has to be passed to the predictor function:\cr
#' \code{x} sample locations (scaled to values between 0 and 1)\cr
#' \code{y} observations at sample locations (see parameters)\cr
#' \code{thetaLower} lower boundary for theta (see parameters)\cr
#' \code{thetaUpper} upper boundary for theta (see parameters)\cr
#' \code{algTheta} algorithm to find theta (see parameters)\cr
#' \code{budgetAlgTheta} budget for the above mentioned algorithm (see parameters)\cr
#' \code{optimizeP} boolean that specifies whether the exponents (\code{p}) were optimized (see parameters)\cr
#' \code{normalizeymin} minimum in normalized space\cr
#' \code{normalizeymax} maximum in normalized space\cr
#' \code{normalizexmin} minimum in input space\cr
#' \code{normalizexmax} maximum in input space\cr
#' \code{dmodeltheta} vector of activity parameters\cr
#' \code{Theta} log_10 vector of activity parameters (i.e. \code{log10(dmodeltheta)})\cr
#' \code{dmodellambda} regularization constant (nugget) \cr
#' \code{Lambda} log_10 of regularization constant (nugget) (i.e. \code{log10(dmodellambda)})\cr
#' \code{yonemu} \code{Ay-ones*mu} \cr
#' \code{ssq} sigma square\cr
#' \code{mu} mean mu\cr
#' \code{Psi} matrix large Psi\cr
#' \code{Psinv} inverse of Psi\cr
#' \code{nevals} number of Likelihood evaluations during MLE
#'
#' @export
#' @seealso \code{\link{predict.cobbsGPR}} 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#'
#' @examples
#' ## Test-function:
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' set.seed(1)
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model with default settings
#' fit <- gaussianProcessRegression(x,y)
#' ## Print model parameters
#' print(fit)
#' ## Predict at new location
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
#' ## 
#' ## Next Example: Handling factor variables
#' ##
#' ## create a test function:
#' braninFunctionFactor <- function (x) {  
#' y <- (x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#' 		10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' if(x[3]==1)
#' 		y <- y +1
#' else if(x[3]==2)
#' 		y <- y -1
#' y  
#' }
#' ## create training data
#' set.seed(1)
#' x <- cbind(runif(50)*15-5,runif(50)*15,sample(1:3,50,replace=TRUE))
#' y <- as.matrix(apply(x,1,braninFunctionFactor))
#' ## fit the model (default: assume all variables are numeric)
#' fitDefault <- gaussianProcessRegression(x,y,control = list(algTheta=DEinterface))
#' ## fit the model (give information about the factor variable)
#' fitFactor <- gaussianProcessRegression(x,y,control = 
#'		list(algTheta=DEinterface,types=c("numeric","numeric","factor")))
#' ## create test data
#' xtest <- cbind(runif(200)*15-5,runif(200)*15,sample(1:3,200,replace=TRUE))
#' ytest <- as.matrix(apply(xtest,1,braninFunctionFactor))
#' ## Predict test data with both models, and compute error
#' ypredDef <- predict(fitDefault,xtest)$y
#' ypredFact <- predict(fitFactor,xtest)$y
#' mean((ypredDef-ytest)^2)
#' mean((ypredFact-ytest)^2)
###################################################################################
gaussianProcessRegression <- function(x, y, control=list()){
  x <- data.matrix(x) #TODO data.matrix is better than as.matrix, because it always converts to numeric, not character! (in case of mixed data types.)
  npar <- length(x[1,])
	con<-list(thetaLower=1e-6, thetaUpper=1e12, 
		types=rep("numeric",npar),
		algTheta=DEinterface, budgetAlgTheta=200, 
		optimizeP= FALSE,
		useLambda=TRUE, lambdaLower = -16, lambdaUpper = 0, 
		startTheta=NULL, reinterpolate=FALSE, target="y")
	con[names(control)] <- control
	control<-con
	
	fit <- control
	k <- ncol(x)
	fit$x <- x
	fit$y <- y
	# normalize input data
	ymin <- 0
	ymax <- 1
	fit$normalizeymin<-ymin
	fit$normalizeymax<-ymax
	res <- normalizeMatrix(fit$x, ymin, ymax)
	fit$scaledx <-res$y
	fit$normalizexmin<-res$xmin
	fit$normalizexmax<-res$xmax
	LowerTheta <- rep(1,k)*log10(fit$thetaLower)

	UpperTheta <- rep(1,k)*log10(fit$thetaUpper)

	#Wrapper for optimizing theta  based on likelihood function
	fitFun <- function (x, fX, fy,optimizeP,useLambda,penval){ #todo vectorize?
		as.numeric(gprLikelihood(x,fX,fy,optimizeP,useLambda,penval)$NegLnLike)
	}
	n<-nrow(fit$x) #number of observations
	
	
	if(is.null(fit$startTheta))
		x1 <-  rep(n/(100*k),k) # start point for theta
	else
		x1 <- fit$startTheta

	A<-matrix(0,k,n*n)

	for(i in 1:k){
		if(control$types[i]!="factor"){
			A[i,]<-as.numeric(as.matrix(dist(fit$scaledx[,i]))) #euclidean distance
		}else	{
			tmp <- outer(fit$scaledx[,i],fit$scaledx[,i],'!=') #hamming distance
			class(tmp) <- "numeric"
			A[i,]<-tmp
		}
	}
	
	if(fit$optimizeP){ # optimize p
		LowerTheta <- c(LowerTheta, rep(1,k)*0.01)
		UpperTheta <- c(UpperTheta, rep(1,k)*2)		
		x3 <- rep(1,k)* 1.9 #start values for p
		x0 <- c(x1,x3)
	}else{ # p  is fixed to 2 and the array A is completely precalculated
		A<-A^2
		x0 <- c(x1)
	}			
	if(fit$useLambda){
		# start value for lambda:
		x2 <-  (fit$lambdaUpper + fit$lambdaLower)/2
		x0 <- c(x0,x2)
		#append regression constant lambda (nugget)
		LowerTheta <- c(LowerTheta,fit$lambdaLower)
		UpperTheta <- c(UpperTheta,fit$lambdaUpper)
	}	
	#force x0 into bounds
	x0 <- pmin(x0,UpperTheta)
	x0 <- pmax(x0,LowerTheta)
  x0 <- matrix(x0,1) #matrix with one row
	opts<-list(funEvals=fit$budgetAlgTheta*length(x0))
	#determine a good penalty value (based on number of samples and variance of y)
	penval <- n*log(var(y)) + 1e4
	res <- fit$algTheta(x=x0,fun=
							function(x,fX,fy,optimizeP,useLambda,penval){
								if(!is.matrix(x)){
									fitFun(x,fX,fy,optimizeP,useLambda,penval)
								}
								else{
									apply(x,1,fitFun,fX,fy,optimizeP,useLambda,penval)
								}
							},
            lower=LowerTheta,upper=UpperTheta, 
						control=opts,fX=A,fy=fit$y,optimizeP=fit$optimizeP,useLambda=fit$useLambda,penval=penval)	
	if(is.null(res$xbest))res$xbest<-x0;
	Params <- res$xbest
	nevals <- as.numeric(res$count[[1]])
	
	fit$dmodeltheta<-10^Params[1:k]
	if(fit$optimizeP){	
		fit$P<-Params[(k+1):(2*k)]
	}
	if(fit$useLambda){
		fit$Lambda <- Params[length(Params)];
		fit$dmodellambda<-10^Params[length(Params)]
	}else{
		fit$Lambda <- -Inf;
		fit$dmodellambda<-0
	}
	# extract model parameters:
	fit$Theta <- Params[1:k]
	res <- gprLikelihood(c(fit$Theta,fit$P, fit$Lambda),A,fit$y,fit$optimizeP,fit$useLambda);
	
	if(is.na(res$Psinv[1]))
		stop("gaussianProcessRegression failed to produce a valid correlation matrix. Consider activating the nugget/regularization constant via useLambda=TRUE in the control list.")
	
	fit$yonemu <- res$yonemu
	fit$ssq <- as.numeric(res$ssq)
	fit$mu <- res$mu
	fit$Psi <- res$Psi
	fit$Psinv <- res$Psinv
	fit$nevals <- nevals
	fit$like <- res$NegLnLike
	fit$returnCrossCor <- FALSE
	
	## calculate observed minimum
	xlist <- split(x, 1:nrow(x))
	uniquex <- unique(xlist)
	ymean <- NULL
	for(xi in uniquex){
		ind <- xlist %in% list(xi)
		ymean <- c(ymean, mean(y[ind]))
	}	
  fit$min <- min(ymean)
		
	class(fit)<- "cobbsGPR"
	fit
}

###################################################################################
#' Normalize design 2
#'
#' Normalize design with given maximum and minimum in input space. Supportive function for GPR model, not to be used directly.
#' 
#' @param x design matrix in input space (n rows for each point, k columns for each parameter)
#' @param ymin minimum vector of normalized space
#' @param ymax maximum vector of normalized space
#' @param xmin minimum vector of input space
#' @param xmax maximum vector of input space
#'
#' @return normalized design matrix
#' @seealso \code{\link{gaussianProcessRegression}}
#' @keywords internal
###################################################################################
normalizeMatrix2 <- function (x,ymin,ymax,xmin,xmax){ 
	rangex <- xmax-xmin
	rangey <- ymax-ymin
	s <- dim(x)[1]
	rangey * (x-matrix(rep(xmin,s),nrow=s,byrow=TRUE))/matrix(rep(rangex,s),nrow=s,byrow=TRUE) + ymin
}

###################################################################################
#' Normalize design
#' 
#' Normalize design by using minimum and maximum of the design values for input space. Supportive function for GPR model, not to be used directly.
#'
#' @param x design matrix in input space
#' @param ymin minimum vector of normalized space
#' @param ymax maximum vector of normalized space
#'
#' @return normalized design matrix
#' @seealso \code{\link{gaussianProcessRegression}}
#' @keywords internal
###################################################################################
normalizeMatrix <- function(x,ymin, ymax){
	# Return the maximum from each row:
	xmax <- apply(x,2,max)
	# Return the minimum from each row:
	xmin <- apply(x,2,min)
	s <- dim(x)[1]
	rangex <- xmax-xmin
	rangey <- ymax-ymin
	xmin[rangex==0] <- xmin[rangex==0]-0.5
	xmax[rangex==0] <- xmax[rangex==0]+0.5
	rangex[rangex==0] <- 1
	y <- rangey * (x-matrix(rep(xmin,s),nrow=s,byrow=TRUE))/matrix(rep(rangex,s),nrow=s,byrow=TRUE) + ymin
	list(y=y,xmin=xmin,xmax=xmax)
}

###################################################################################
#' Calculate negative log-likelihood
#' 
#' Used to determine theta/lambda values for the GPR model in \code{\link{gaussianProcessRegression}}.
#' Supportive function for GPR model, not to be used directly.
#'
#' @param x vector, containing parameters log10(theta), log10(lambda) and p.
#' @param AX 3 dimensional array, constructed by gaussianProcessRegression from the sample locations
#' @param Ay vector of observations at sample locations
#' @param optimizeP boolean, whether to optimize parameter p (exponents) or fix at two.
#' @param useLambda boolean, whether to use nugget
#' @param penval a penalty value which affects the value returned for invalid correlation matrices / configurations
#'
#' @return list with elements\cr
#' \code{NegLnLike}  concentrated log-likelihood *-1 for minimising \cr
#' \code{Psi} correlation matrix\cr
#' \code{Psinv} inverse of correlation matrix (to save computation time during prediction)\cr
#' \code{mu} \cr
#' \code{ssq}
#' @seealso \code{\link{gaussianProcessRegression}}
#' @export
#' @keywords internal
###################################################################################
gprLikelihood <- function(x,AX,Ay,optimizeP=FALSE,useLambda=TRUE,penval=1e8){
	## todo: this penalty value should not be a hard constant.
	## the scale of the likelihood (n*log(SigmaSqr) + LnDetPsi)
	## at least depends on log(var) and the number of samples
	## hence, large number of samples may lead to cases where the
	## penality is lower than the likelihood of most valid parameterizations
	## suggested solution (still requires testing): 
	#penval <- n*log(var(y)) + 1e4
	
	
	nx<-nrow(AX)
	theta <- 10^x[1:nx]
	if(optimizeP){	
		AX <- abs(AX)^(x[(nx+1):(2*nx)])
	}
	lambda <- 0
	if(useLambda){
		lambda <- 10^x[length(x)]
	}
	if( any(theta==0) ||  any(is.infinite(c(theta,lambda)))){ #for instance caused by bound violation
		penalty <- penval 
		return(list(NegLnLike=penalty,Psi=NA,Psinv=NA,mu=NA,ssq=NA))
	}
	n <- dim(Ay)[1]
	
	Psi <- exp(-matrix(colSums(theta*AX),n,n))
	if(useLambda){
		Psi <- Psi+diag(lambda,n)
	}
	
	if(any(is.infinite(Psi))){ # this is required especially if distance matrices are forced to be CNSD/NSD and hence have zero distances
	  penalty <- penval 
		return(list(NegLnLike=penalty,Psi=NA,Psinv=NA,mu=NA,SSQ=NA)) #todo: add something like smallest eigenvalue to penalty?
	}	
		
	## Check whether Psi is ill-conditioned
  kap <- rcond(Psi)  
	if(is.na(kap)){
		kap <- 0
	}
	if(kap < 1e-10){ 
		#warning("Correlation matrix is ill-conditioned (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty.")
		return(list(NegLnLike=penval-log10(kap),Psi=NA,Psinv=NA,mu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=TRUE))
	}	
		
	## cholesky decomposition
	cholPsi <- try(chol(Psi), TRUE) 

	## give penalty if fail
	if(class(cholPsi)[1] == "try-error"){
		#warning("Correlation matrix is not positive semi-definite (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty.")
		penalty <- penval - log10(min(eigen(Psi,symmetric=TRUE,only.values=TRUE)$values))
		return(list(NegLnLike=penalty,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
		
	#calculate natural log of the determinant of Psi (numerically more reliable and also faster than using det or determinant)
	LnDetPsi <- 2*sum(log(abs(diag(cholPsi))))
	
	#inverse with cholesky decomposed Psi
	Psinv <- try(chol2inv(cholPsi), TRUE) 
	
  ## give penalty if failed
	if(class(Psinv)[1] == "try-error"){
		#warning("Correlation matrix is not positive semi-definite (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty.")
		penalty <- penval - log10(min(eigen(Psi,symmetric=TRUE,only.values=TRUE)$values))
		return(list(NegLnLike=penalty,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
  
	psisum <- sum(Psinv) #this sum of all matrix elements may sometimes become zero, which may be caused by inaccuracies. then, the following may help
	if(psisum==0){
		psisum <- as.numeric(rep(1,n) %*% Psinv %*% rep(1,n))
		if(psisum==0){ #if it is still zero, return penalty
			#warning("Sum of elements in inverse correlation matrix is zero (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty.")
			return(list(NegLnLike=penval,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
		}
	}		
	
	mu <- sum(Psinv%*%Ay)/psisum
	if(is.infinite(mu)|is.na(mu)){
		#warning("MLE estimate of mu is infinite or NaN (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty.")
		return(list(NegLnLike=penval,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}			
	
	yonemu <- Ay-mu 
	SigmaSqr <- (t(yonemu)%*%Psinv%*%yonemu)/n
	if(SigmaSqr < 0){
		#warning("Maximum Likelihood Estimate of model parameter sigma^2 is negative (During Maximum Likelihood Estimation in gaussianProcessRegression). Returning penalty. ")
		return(list(NegLnLike=penval-SigmaSqr,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
	
	NegLnLike <- n*log(SigmaSqr) + LnDetPsi
  if(is.na(NegLnLike)|is.infinite(NegLnLike)){
    return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}
	list(NegLnLike=NegLnLike,Psi=Psi,Psinv=Psinv,mu=mu,yonemu=yonemu,ssq=SigmaSqr)
}

###################################################################################
#' Predict GPR Model
#' 
#' Predict with GPR model produced by \code{\link{gaussianProcessRegression}}.
#'
#' @param object GPR model (settings and parameters) of class \code{cobbsGPR}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted mean \code{y}, uncertainty / standard deviation \code{s} (optional) and expected improvement \code{ei} (optional). 
#' Whether \code{s} and \code{ei} are returned is specified by the vector of strings \code{object$target}, which then contains \code{"s"} and \code{"ei"}.
#'
#'
#' @examples
#' ## Test-function:
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 +
#'  10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model
#' fit <- gaussianProcessRegression(x,y)
#' fit$target <- c("y","s","ei")
#' ## first estimate error with regressive predictor
#' predict(fit,x)
#'
#' @seealso \code{\link{gaussianProcessRegression}}
#' @export
#' @keywords internal
###################################################################################
predict.cobbsGPR <- function(object,newdata,...){
	x<-newdata
	if(object$reinterpolate){
		return(predictGPRReinterpolation(object,x,...))
	}
  
	x <- repairNonNumeric(x,object$types) #round to nearest integer.
	
  #normalize input x
	x <- normalizeMatrix2(data.matrix(x),0,1,object$normalizexmin,object$normalizexmax)
	AX <- object$scaledx
	theta <- object$dmodeltheta
	Psinv <- object$Psinv #fixed: does not need to be computed, is already done in likelihood function
	n <- dim(AX)[1]
	mu <- object$mu
	yonemu <- object$yonemu	
	SigmaSqr <- object$ssq
 
  k <- nrow(x)
	nvar <- ncol(x)
	
	if(object$optimizeP)
	  p <- object$P
	else
	  p <- rep(2,nvar)	
		
  psi <- matrix(0,k,n)
	for (i in 1:nvar){ #todo nnn number variables 
	  #psi[,i] <- colSums(theta*(abs(AX[i,]-t(x))^p))
	  tmp <- expand.grid(AX[,i],x[,i])
    if(object$types[i]=="factor"){
      tmp <- as.numeric(tmp[,1]!=tmp[,2])^p[i]
    }else{
      tmp <- abs(tmp[,1]-tmp[,2])^p[i]
	  }
    psi <- psi + theta[i] * matrix(tmp,k,n,byrow = T)	  
	}
  
	psi <- exp(-psi)
	f <- as.numeric(psi%*%(Psinv%*%yonemu))+mu #vectorised			#TODO: Psinv %*% yonemu can be computed ahead of time
	##########################################################################
	res <- list(y=f)
	if (any(object$target %in% c("s","ei"))){
		lambda <- object$dmodellambda;
		SSqr <- SigmaSqr*(1+lambda-diag(psi%*%(Psinv%*%t(psi)))) 
		s <- sqrt(abs(SSqr))
    res$s <- s
    if(any(object$target == "ei")){
      res$ei <- expectedImprovement(f,s,object$min)
    }    
	}
	if(object$returnCrossCor)
		res$psi <- psi
  res
}


###################################################################################
#' Predict GPR Model (Re-interpolating)
#' 
#' GPR predictor with re-interpolation to avoid stalling the optimization process which employs this model as a surrogate.
#' This is supposed to be used with deterministic experiments, which do need a non-interpolating model that avoids predicting non-zero error at sample locations.
#' This can be useful when the model is deterministic (i.e. repeated evaluations of one parameter vector do not yield different values) but does have a "noisy" structure (e.g. due to computational inaccuracies, systematical error).
#'
#' Please note that this re-interpolation implementation will not necessarily yield values of exactly zero at the sample locations used for model building. Slight deviations can occur.
#'
#' @param object GPR model (settings and parameters) of class \code{cobbsGPR}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted mean \code{y}, uncertainty \code{s} (optional) and expected improvement \code{ei} (optional). 
#' Whether \code{s} and \code{ei} are returned is specified by the vector of strings \code{object$target}, which then contains "s" and "ei.
#'
#' @examples
#' ## Test-function:
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model
#' fit <- gaussianProcessRegression(x,y,control=list(reinterpolate=FALSE))
#' fit$target <- c("y","s")
#' ## first estimate error with regressive predictor
#' sreg <- predict(fit,x)$s
#' ## now estimate error with re-interpolating predictor
#' sreint <- predictGPRReinterpolation(fit,x)$s
#' ## equivalent:
#' fit$reinterpolate <- TRUE
#' sreint2 <- predict(fit,x)$s
#' print(sreg)
#' print(sreint)
#' print(sreint2)
#' ## sreint should be close to zero, significantly smaller than sreg
#'
#' @seealso \code{\link{gaussianProcessRegression}}, \code{\link{predict.cobbsGPR}}
#' @export
#' @keywords internal
###################################################################################
predictGPRReinterpolation <- function(object,newdata,...){
	x <- newdata
	x <- repairNonNumeric(x,object$types) #round to nearest integer.
	#normalize input x
	x <- normalizeMatrix2(data.matrix(x),0,1,object$normalizexmin,object$normalizexmax)
	AX <- object$scaledx
	theta <- object$dmodeltheta
	lambda <- object$dmodellambda
	Psi <- object$Psi
	Psinv <- object$Psinv #fixed: does not need to be computed, is already done in likelihood function
	n <- dim(AX)[1]
	mu <- object$mu
	yonemu <- object$yonemu	
	#
	PsiB <- Psi-diag(lambda,n)+diag(.Machine$double.eps,n) #TODO this and the following can be precomputed during model building
	SigmaSqr <- as.numeric(t(yonemu)%*%Psinv%*%PsiB%*%Psinv%*%yonemu)/n
	#	
	k <- nrow(x)
	nvar <- ncol(x)
	
	if(object$optimizeP)
	  p <- object$P
	else
	  p <- rep(2,nvar)	
	
	psi <- matrix(0,k,n)
	for (i in 1:nvar){ #todo nnn number variables 
	  #psi[,i] <- colSums(theta*(abs(AX[i,]-t(x))^p))
	  tmp <- expand.grid(AX[,i],x[,i])
	  if(object$types[i]=="factor"){
	    tmp <- as.numeric(tmp[,1]!=tmp[,2])^p[i]
	  }else{
	    tmp <- abs(tmp[,1]-tmp[,2])^p[i]
	  }
	  psi <- psi + theta[i] * matrix(tmp,k,n,byrow = T)	  
	}
	
	
	psi <- exp(-psi)
	f <- as.numeric(psi%*%(Psinv%*%yonemu))+mu #vectorised
	##########################################################################
	res <- list(y=f)
	if (any(object$target %in% c("s","ei"))){
		Psinv <- try(solve.default(PsiB), TRUE) 
		if(class(Psinv)[1] == "try-error"){
			Psinv<-ginv(PsiB)
		}	
		#
		SSqr <- SigmaSqr*(1-diag(psi%*%(Psinv%*%t(psi)))) #vectorised
		s <- sqrt(abs(SSqr))
    res$s <- s
    if(any(object$target == "ei")){
      res$ei <- expectedImprovement(f,s,object$min)
    }    
	}
	if(object$returnCrossCor)
		res$psi <- psi
  res
}



###################################################################################
#' Print Function GPR
#'
#' Print information about a GPR model fit, as produced by \code{\link{gaussianProcessRegression}}.
#'
#' @rdname print
#' @method print cobbsGPR
#  @S3method print cobbsGPR
#' @param x	fit returned by \code{\link{gaussianProcessRegression}}.
#' @param ... additional parameters	
#' @export
#' @keywords internal
###################################################################################
print.cobbsGPR <- function(x,...){
	cat("------------------------\n")
	cat("GPR model.\n")
	cat("------------------------\n")
	cat("Estimated activity parameters (theta) sorted \n")
	cat("from most to least important variable \n")
	cat(paste("x",order(x$dmodeltheta,decreasing=TRUE),sep="",collaps=" "))
	cat("\n")	
	cat(sort(x$dmodeltheta,decreasing=TRUE))	
	cat("\n \n")
	cat("exponent(s) p:\n")
	if(x$optimizeP)
		cat(x$P)
	else
		cat(2)
	cat("\n \n")
	cat("Estimated regularization constant (or nugget) lambda:\n")	
	cat(x$dmodellambda)
	cat("\n \n")
	cat("Number of Likelihood evaluations during MLE:\n")	
	cat(x$nevals)
	cat("\n")	
	cat("------------------------\n")
}

