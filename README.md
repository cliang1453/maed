<h1 align="center">mAED</h1>
<h4 align="center">R Library for Multi-Stage Adaptive Enrichment Design</h4>

___mAED___ (Multi-Stage Adaptive Enrichment Design) implements a solver with user API for nonsmooth dual problems under highly sparse and structured settings in 2AED and mAED problems. The core algorithm is implemented in C++ with Eigen3, a portable library for high performance linear algebra, and with package [PRIMAL](https://github.com/ShenQianli/primal)[1], a unified framework of parametric simplex methods for sparse learning. 

## Table of contents

- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Directory structure](#directory-structure)
- [Installation](#installation)
- [Performance](#performance)
- [References](#references)

## Introduction

The current approaches formulate mAED problem as a general linear programming (LP) problem [2], which is computationally expensive to solve. This library is proposed for new customized algorithms with three key features: 1) It provides a highly efficient solver to tackle a large and important class of LP problems; 2) It provides a solution for multi-stage decision-making problems with Bayes risk constraints; 3) It provides additional functions such as visualizing the optimal decision maps. 

# We also provide tutorials on the theoretical background and the code. Please see ``tutorial`` folder for tutorials.

## Directory structure
# The directory is organized as follows:
# * [__src__](src): C++ implementation of the PSM algorithm.
# 	* [__api.cpp__](api.cpp): C API as an interface for R and Python package.
# 	* [__PSM.cpp__](PSM.cpp): Core implemetation of the PSM solver. 
# * [__include__](include)
# 	* [__PSM__](PSM): declarations of the C++ implementation
# 	* [__Eigen__](eigen3): Eigen3 header files for high performance linear algebra.
# * [__R-package__](R-package): R wrapper for the source code.
# * [__python-package__](python-package): Python wrapper for the source code.
# * [__tutorials__](tutorials): tutorials for using the code in R and Python.
# * [__profiling__](profiling): profiling the performance from R package.
# * [__Makefile__](Makefile):Makefile local configurations.
# * [__CmakeLists.txt__](CmakeLists.txt):Makefile local configurations.


## Installation
# Third-party dependencies. The installation only depends on Eigen3's header files which are included in this github repo as a submodule. We use Eigen3 as an independent fully portable module so any existing Eigen3 installation will not have conflict with PSM installation. There are two ways to install the mAED R package.
# - Installing from CRAN (recommended). The R package is hosted on CRAN. The easiest way to install R package is by running the following command in R
# ```R
# install.packages("maed")
# ```
# - Installing from source code.
# ```bash
# $ git clone --recurse-submodules https://github.com/cliang1453/maed.git
# $ cd maed; make Rinstall
# ```

## Examples

# Now we illustrate the R user interface using Dantzig selector as an example.
#  ```R
#  set.seed(1024)
#  library(PRIMAL)
#  ## Generate the design matrix and coefficient vector
#  n <- 100; d <- 250; c <- 0.5
#  # n sample number, d dimension, c correlation parameter
#  X <- scale(matrix(rnorm(n*d),n,d)+c*rnorm(n))/sqrt(n-1)*sqrt(n)
#  s <- 5 # sparsity level
#  beta <- c(runif(s,-1,1), rep(0, d-s))
#  Y <- X%*%beta + rnorm(n)
#  ## Dantzig selection solved with parametric simplex method
#  fit.dantzig <- Dantzig_solver(X, Y, max_it = 100, lambda_threshold = 0.01)
#  ## print lambdas used and number of nonzero coefficients for each lambda
#  print(fit.dantzig$lambda)
#  print(fit.dantzig$df)
#  ## Visualize the solution path
#  plot.primal(fit.dantzig)
# ```

## Performance
# ```bash
# $cd profiling
# $Rscript benchmark.R
# ```
# 
# We compare the timing performance of our package with R package "fastclime". We fix the sample size n to be 200 and vary the data dimension d from 1000 to 7000 and 200 to 1000. Each entries of X is independent Gaussian and Gaussianized such that the column has uniform norm. We randomly select 2% entries from vector θ to be nonzero. Algorithm will stop when λ is less than <img src="http://chart.googleapis.com/chart?cht=tx&chl= $$2\sqrt{log(d)/n}$$" style="border:none;">. 
# - Dantzig selector. PRIMAL achieves similar optimization performance to fastclime. But PRIMAL is 6 to 10 times faster than fastclime.
# - Compressed sensing. PRIMAL is 2 times faster than fastclime and achieves similar optimization.
# ![Performance_R](https://github.com/ShenQianli/primal/blob/master/profiling/images/performance_R.jpg)

## References

[1] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.

[2] M. Rosenblum, E. X. Fang, and H. Liu, Optimal, two stage, adaptive enrichment designs for randomized trials, using sparse linear programming, Under minor revision at Journal of Royal Statistical Society, Series B, 2017.


