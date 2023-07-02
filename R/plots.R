###################################################################################
#' Validate Benchmark plot
#' 
#' This plot is intended to plot the error (over runtime) for a benchmark function
#' used to test optimization algorithms. For details on usage, please see the example.
#' 
#' 
#' @param groundtruth a data frame created with \code{\link{loggedExperiment}}. This data set represents the groundruth, which means that the experiment was performed with the 'true' objective function.
#' @param data a list of data frames created with \code{\link{loggedExperiment}}. Each data frame represents the result of an experiment with a different emulation of the ground-truth function.
#' @param names a character vector. Each represents a name for the method corresponding to an element of the \code{data} list.
#'
#' @return a ggplot instance
#' 
#' @examples
#' ## generate some example data, with groundtruth sqrt(|x|), and emulations x^2, log(|x|+0.1).
#' lower <- -1
#' upper <- 1
#' require(nloptr)
#' expr <- expression(
#'   res <- lbfgs(x0=runif(1,lower,upper),fn = fnlog,
#'                      lower=lower,upper=upper)
#' )
#' res1 <- loggedExperiment(expr, function(x)sqrt(abs(x)), 1:10)
#' res2 <- loggedExperiment(expr, function(x)x^2, 1:10)
#' res3 <- loggedExperiment(expr, function(x)log(abs(x)+0.1), 1:10)
#' ## plot 
#' p <- plotBenchmarkValidation(res1,list(res2,res3),names=c("square","log"))
#' print(p)
#'
#' @export
###################################################################################
plotBenchmarkValidation <- function(groundtruth, data, names=NULL){
  if(is.null(names))
    names <- letters[1:length(data)]
  df <- groundtruth[,c("it","seed","cy")]
  for(i in 1:(length(data))){
    #merge
    df <- merge(df,data[[i]][,c("it","seed","cy")],by=c("it","seed"),all=TRUE)
  }
  df <- df[order(df$seed),]
  df <- data.frame(nafill(df,"locf"))
  names(df) <- c("it","seed","groundtruth",names)
  ## scale
  df[,-(1:2)] <- apply(df[,-(1:2)],2,scale01)
  ## get absolute deviation
  nmserr <- paste(names,"err",sep=".")
  ddf <- NULL
  for(i in 1:(length(data))){
    #merge #[[nmserr[i]]]
    dftmp <- data.frame(it=df$it)
    dftmp[[nmserr[i]]] <-  abs(df$groundtruth - df[[names[i]]])
    dftmp <- dftmp %>% group_by(.data$it) %>%  
      summarise(upper = quantile(!!sym(nmserr[i]),probs=0.75),
                lower = quantile(!!sym(nmserr[i]),probs=0.25),
                med = median(!!sym(nmserr[i]))) 
    dftmp$method=names[i]
    ddf <- rbind(ddf,dftmp)
  }
  p <- ggplot(
    ddf,
    aes_string(x="it",y="med",group="method"))  + 
    geom_ribbon(aes_string(ymin="lower", ymax="upper",fill="method"), alpha=0.5) +
    geom_line(aes_string(color="method"))+ theme_classic() +
    ylab("error")+xlab("evaluations")+
    theme(
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank(),
      strip.background = element_blank()
    )+ theme(legend.position="top")
  return(p)
}

###################################################################################
#' Performance plot
#' 
#' Create a simple performance plot, showing the scaled performance of an algorithm over time (based on experiments with different objective functions).
#' 
#' @param data a list of data frames created with \code{\link{loggedExperiment}}. Each data frame represents the result of an experiment with a different objective function.
#' @param names a character vector. Each represents a name for the underlying objective function of the \code{data} list. 
#'
#' @return a ggplot instance
#' 
#' @examples
#' ## create data by running optim L-BFGS-B with two objective functions
#' lower <- -1
#' upper <- 1
#' require(nloptr)
#' expr <- expression(
#'   res <- lbfgs(x0=runif(1,lower,upper),fn = fnlog,
#'                      lower=lower,upper=upper)
#' )
#' res1 <- loggedExperiment(expr, function(x)x^2, 1:10)
#' res2 <- loggedExperiment(expr, function(x)sqrt(abs(x)), 1:10)
#' ## create the plot from the results
#' p <- plotBenchmarkPerformance(list(res1,res2),names=c("square","root"))
#' print(p)
#'
#' @export
###################################################################################
#res1 <- rbind(data.frame(it=c(1:10),seed=1,y=NA,cy=10:1),
#              data.frame(it=c(1:10),seed=2,y=NA,cy=(10:1)^1.1),
#              data.frame(it=c(1:10),seed=3,y=NA,cy=(10:1)^1.2)
#              )
#res2 <- rbind(data.frame(it=c(1:10),seed=1,y=NA,cy=(10:1)^2),
#              data.frame(it=c(1:10),seed=2,y=NA,cy=(10:1)^2.1),
#              data.frame(it=c(1:10),seed=3,y=NA,cy=(10:1)^1.8)
#              )
#p <- plotBenchmarkPerformance(list(res1,res2))
#print(p)
plotBenchmarkPerformance <- function(data,names=NULL){
  if(is.null(names))
    names <- letters[1:length(data)]
  df <- data[[1]][,c("it","seed","cy")]
  for(i in 2:(length(data))){
    #merge
    df <- merge(df,data[[i]][,c("it","seed","cy")],by=c("it","seed"),all=TRUE)
  }
  df <- df[order(df$seed),]
  df <- data.frame(data.table::nafill(df,"locf"))
  names(df) <- c("it","seed",names)
  ## scale
  df[,-(1:2)] <- apply(df[,-(1:2)],2,scale01)
  ddf <- reshape2::melt(df,id.vars=c("it","seed"))
  names(ddf) <- c("it","seed","method","value")
  ddf <- ddf %>% group_by(.data$it,.data$method) %>%  
    summarise(upper = quantile(.data$value,probs=0.75),
              lower = quantile(.data$value,probs=0.25),
              med = median(.data$value)) 
  p <- ggplot(
    ddf,
    aes_string(x="it",y="med",group="method"))  + 
    geom_ribbon(aes_string(ymin="lower", ymax="upper",fill="method"), alpha=0.5) +
    geom_line(aes_string(color="method"))+ theme_classic() +
    ylab("y")+xlab("evaluations")+
    theme(
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank(),
      strip.background = element_blank()
    )+ theme(legend.position="top")
  return(p)
}