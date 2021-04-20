# COBBS: Continuous Optimization Benchmarks By Simulation 

## Brief Description
This is a package for R, providing some tools to create simulation-based benchmarks for continuous
optimization, and to run simple experiments.

## Installation

From your R console run

`devtools::install_github("martinzaefferer/COBBS")`

Note, that this requires the devtools package being installed (`install.packages("devtools")`).

## Use

See the example in

`require(COBBS)` 

`?generateCOBBS` 

## Simple example

```
require(COBBS)

# objective  function:
fun <- function(x){
  sin(33*x) + sin(49*x - 0.5)+ x 
}
# "vectorize" target
groundtruth <- function(x){sapply(x,fun)}


# samples (training data)
x <- cbind(c(0.13,0.6,0.62,0.67,.75,0.79,0.8,0.86,0.9,0.95,0.98))
# corresponding training observations
y <- cbind(groundtruth(x))

## specify some model configuration
mc  <- list(useLambda=FALSE,thetaLower=1e-6,thetaUpper=1e12)
## and some configuration details for the simulation
cntrl <- list(modelControl=mc,
              nsim=1,
              seed=1,
              method="spectral",
              Ncos = 10, #for simplicity sake a relatively low number of cosines
              #method="decompose",xsim=matrix(c(runif(200,-32,32)),,1),
              conditionalSimulation=TRUE
)
## generate model and functions
cobbsResult <- generateCOBBS(x,y,cntrl)
xplot <- seq(from=0,to=1, length.out=100)
## plot ground truth function
plot(xplot,groundtruth(xplot),type="l",ylim=c(-2,4),lty=2)
## plot estimation (prediction) 
lines(xplot,cobbsResult$estimation(xplot),col="red")
## plot simulation
lines(xplot,cobbsResult$simulation[[1]](xplot),col="blue")
```

Now you might be interested in what the "simulated" function (blue line) actually looks like. 
It's general shape is a sum of cosines, with parameters omega and phi:
``` 
environment(cobbsResult$simulation[[1]])$fun
```
Note, that the function first performs a normalization step (so that input x scales from 0 to 1).
If desired, you can extract the parameters from the environment:
``` 
environment(environment(cobbsResult$simulation[[1]])$fun)$omega
environment(environment(cobbsResult$simulation[[1]])$fun)$phi
```
Both of these parameters are produced stochastically, but depend on the trained model.
In addition, a multiplier (depends on Gaussian process' global variance and number of cosines Ncos)
and the mean (of the Gaussian process) are used.

``` 
environment(environment(cobbsResult$simulation[[1]])$fun)$multiplier
environment(environment(cobbsResult$simulation[[1]])$fun)$mu
```

## Demo
The following is a more complete demonstration, for a 2-dimensional function,
including actual benchmark runs with an optimization algorithm.

```
require(COBBS)
## generate some data
require(smoof)
seed <- 1234
fnbbob <- makeBBOBFunction(dimensions = 2, fid = 23, iid = 2)
groundtruth <- function(x){
  x=matrix(x,,2) 
  apply(x,1,fnbbob)
}
lower = getLowerBoxConstraints(fnbbob)
upper = getUpperBoxConstraints(fnbbob)
dimension <- length(lower)
set.seed(seed)

## prepare an expression that will be run during the experiments
## here: DE
expr <- expression(
  res <- DEinterface(fun = fnlog,lower=lower,upper=upper,control=list(funEvals=dimension*100,populationSize=dimension*20))
)
## run an experiments, with logging
require(COBBS)
resgt <- loggedExperiment(expr, groundtruth, 1,logx = TRUE)
resgt <- resgt[1:(dimension*50),]
x <- as.matrix(resgt[,c(4,5)])
y <- as.matrix(resgt[,2,drop=F])

## specify some model configuration
mc  <- list(useLambda=FALSE,thetaLower=1e-6,thetaUpper=1e12)
## and some configuration details for the simulation
cntrl <- list(modelControl=mc,
              nsim=1,
              seed=seed,
              method="spectral",
              Ncos = 100*dimension,
              #method="decompose",xsim=matrix(c(runif(200,-32,32)),,1),
              conditionalSimulation=TRUE
)
## generate model and functions
cobbsResult <- generateCOBBS(x,y,cntrl)
cobbsResult$fit

## plot trained model (predictor/estimatation and simulation), using SPOT package
SPOT:::plotFunction(groundtruth,lower,upper)
SPOT:::plotFunction(cobbsResult$estimation,lower,upper)
SPOT:::plotFunction(cobbsResult$simulation[[1]],lower,upper)

## prepare an expression that will be run during the experiments
## here: DE
require(nloptr)
expr <- expression(
  res <- DEinterface(fun = fnlog,lower=lower,upper=upper,control=list(funEvals=1000*dimension))
)
## run the experiments, with logging
## with each objective function produced by COBBS

resgt <- loggedExperiment(expr, groundtruth, 1:10,10)
reses <- loggedExperiment(expr, cobbsResult$estimation, 1:10,10)
ressi <- loggedExperiment(expr, cobbsResult$simulation[[1]], 1:10,10)
## plot results
print(plotBenchmarkPerformance(list(resgt,reses,ressi),c("groundtruth","estimation","simulation")))
## plot error, comparing against groundtruth
print(plotBenchmarkValidation(resgt,list(reses,ressi),c("estimation","simulation")))
```
Note, that for simplicity's sake this is not quite the same experiment as in the paper.
To recreate that experiment would require running all three tested algorithms (Nelder-Mead and random search in addition to DE).
Also, this is an experiment repeated on the same problem instances with different initial seeds for each algorithm run.
Wheras in the paper, only a single run is performed for each instance (but experiments are instead repeated by
testing with different BBOB instances (for the respective BBOB function)).

