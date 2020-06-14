################################################################################
#   Copyright (c) 2020 by Martin Zaefferer, Cologne University of Applied Sciences (TH Koeln)
################################################################################
##	This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

################################################################################
#' Continuous Optimization Benchmarks By Simulation
#'
#' This package implements and demonstrates
#' the use of simulation with Gaussian process regression models
#' for benchmarking of continuous optimization algorithms.
#'
#' \tabular{ll}{
#' Package: \tab COBBS\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.0\cr
#' Date: \tab 2020-04-23\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name COBBS-package
#' @aliases COBBS
#' @docType package
#' @title Continuous Optimization Benchmarks By Simulation
#' @author Martin Zaefferer \email{mzaefferer@@gmail.com} 
# @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @keywords package
#' @seealso Interface of main function: \code{\link{generateCOBBS}} 
#' @import DEoptim
#' @import reshape2
#' @import smoof
#' @import ggplot2
#' @import dplyr
#' @import nloptr
#' @importFrom data.table nafill
#' @importFrom stats median quantile
#' @importFrom utils tail
#
NA #ends description
################################################################################