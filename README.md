<h1 align="center">mAED</h1>
<h4 align="center">R Library for Multi-Stage Adaptive Enrichment Design</h4>

___mAED___ (Multi-Stage Adaptive Enrichment Design) implements a solver interface for nonsmooth dual problems under highly sparse and structured settings in mAED problems. The core implementation is based on Eigen3 library[3], a portable C++ library for high performance linear algebra, and [PRIMAL](https://github.com/ShenQianli/primal)[1] package, a unified framework of parametric simplex methods for sparse learning. 

## Table of contents

- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Installation](#installation)
- [Examples](#examples)
- [References](#references)

## Introduction

Multi-stage Adaptive Enrichment Design (mAED) problem is mostly formulated as a general linear programming (LP) problem [2], which is computationally expensive to solve. This library is proposed for implementation of new customized algorithms with several key features: 1) It provides a highly efficient solver to tackle a large and important class of LP problems; 2) It provides a solution for multi-stage decision-making problems with Bayes risk constraints. This package currently supports solving the 2-stage 2-subpopulation AED problem (2AED), with the core optimization engine supported by [PRIMAL](https://github.com/ShenQianli/primal)[1]. This package provides an interface to the LP solver, where the users could specify the 2AED problem settings, e.g., prior distributions, stage settings, subpopulation settings and additional constraints. 

## Installation

### Installing from GitHub

First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type

```
install.packages(devtools)
```

Then load the devtools package and install maed

```
library(devtools)
install_github("cliang1453/maed")
library(maed)
```

## Examples

```R
library(maed) 
n_stage <- 2
n_ppl <- 2
prop_ppl <- c(0.3 0.7)
R_sizes <- c(100, 100)
G_size <- 500
lambda <- "normal"
solver <- "primal"

# User specify the values of D, G, G_c2, loss_table, alpha, beta, J.
est = maed(n_stage=n_stage, n_ppl=n_ppl, prop_ppl=prop_ppl, 
            R_sizes=R_sizes, G_size=G_size,
            D=D, G=G, G_c2=G_c2, alpha=alpha, J=J, beta=beta,
            lambda_0=lambda, loss_tables_0=loss_table,
            lambda_c2=lambda, loss_tables_c2=loss_table,
            solver="primal")

# Visualize the solution path  
plot(est)
```
## What's next
- This package will be in future extended to support solving the general Multi-stage AED (mAED) problem with more advanced methods including row generation and subspace optimization[2], and will be released through CRAN.
- This package will be in future support the visualization of the optimal decision map.

## References

[1] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.  
[2] M. Rosenblum, E. X. Fang, and H. Liu, Optimal, two stage, adaptive enrichment designs for randomized trials, using sparse linear programming, Under minor revision at Journal of Royal Statistical Society, Series B, 2017.  
[3] D. Eddelbuettel and R. François, “Rcpp: Seamless R and C++ integration,” Journal of Statistical Software, vol. 40, no. 8, pp. 1–18, 2011.  
[4] Guennebaud, G., Jacob, B. et al.(2010).  Eigen v3.  http://eigen.tuxfamily.org.  
