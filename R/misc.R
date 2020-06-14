
###################################################################################################
#' Model Fit to Function
#'
#' This function produces an objective function with y=f(x) from a provided model fit.
#' Important note: this function expects \code{predict(object,newdata)} to return
#' a list. The \code{object$target} parameter is a string that determins which list item
#' is returned by the created function. If not set (NULL), \code{object$target} is set to \code{"y"}.
#'
#' @param object fit created by a modeling function, e.g., \code{\link{gaussianProcessRegression}}
#' @param infillCriterion optional parameter, a function that accepts prediction results and a model object. The function should use these
#' to alter the prediction result in a user desired way. For example turning the prediction results of a kriging model (mean and sd) into the expected 
#' improvement criterion
#' @return a function in the style of \code{y=f(x)}, which uses the fitted object to predict \code{y} for sample \code{x}. 
#' @export
#' @keywords internal
###################################################################################################
fit2function <- function(object, infillCriterion = NULL){
    if(is.null(object$target)) 
        object$target <- "y"
    
    evalModelFun <- function(x){  
        res <- predict(object=object,newdata=x)[object$target]
        if(length(res) == 1){
            res <- res[[1]]
        }
        return(res)
    }
    
    if(is.null(infillCriterion)){
        return(
            #evalModelFun
            function(x){
                res <- evalModelFun(x)
                if(is.list(res)){
                    return(res[[1]])
                }
                return(res)
            }
        )
    }
    
    return(
        function(x){
            return(infillCriterion(evalModelFun(x), object))
        }
    )
}

###################################################################################################
#' Repair Non-numeric Values
#'
#' Round non-numeric columns of a matrix, specified by a vector of data given data types.
#'
#' @param x matrix to be rounded
#' @param types data types of the respective columns, numeric columns are specified by \code{"numeric"}.
#' @export
#' @keywords internal
#' @examples
#' x <- matrix(10*runif(12),4,3)
#' types <- c("numeric","factor","factor")
#' repairNonNumeric(x,types)
###################################################################################################
repairNonNumeric <- function(x,types){
	for (i in 1:length(types)){
    if(types[i] != "numeric") #use rounding if not numeric. note that categorical parameters are mapped to integers.
			x[,i] <- round(x[,i])
	}
	x
}

###################################################################################################
#' Expected Improvement
#'
#' Compute the negative logarithm of the Expected Improvement of a set of candidate solutions.
#' Based on mean and standard deviation of a candidate solution,
#' this estimates the expectation of improvement. Improvement
#' considers the amount by which the best known value (best observed value)
#' is exceeded by the candidates.
#'
#' @param mean vector of predicted means of the candidate solutions.
#' @param sd vector of estimated uncertainties / standard deviations of the candidate solutions.
#' @param min minimal observed value.
#' 
#' @return a vector with the negative logarithm of the expected improvement values, -log10(EI).
#'
#' @export
#' @examples
#' mean <- 1:10 #mean of the candidates
#' sd <- 10:1 #st. deviation of the candidates
#' min <- 5 #best known value
#' EI <- expectedImprovement(mean,sd,min)
#' EI
###################################################################################################
expectedImprovement <- function(mean,sd,min){ #NegLogExpImp 
	EITermOne=(min-mean)*pnorm((min-mean)/sd)
	EITermTwo=sd*(1/sqrt(2*pi))*exp(-(1/2)*((min-mean)^2/(sd^2)))
	-log10(EITermOne+EITermTwo+(.Machine$double.xmin))
}