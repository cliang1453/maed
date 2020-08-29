<h1 align="center">mAED</h1>
<h4 align="center">R Library for Multi-Stage Adaptive Enrichment Design</h4>

___mAED___ (Multi-Stage Adaptive Enrichment Design) implements a solver for nonsmooth dual problems under highly sparse and structured settings in 2AED problems. The core algorithm is implemented in C++ with Eigen3, a portable library for high performance linear algebra, and with package [PRIMAL](https://github.com/ShenQianli/primal)[1], a unified framework of parametric simplex methods for sparse learning. 

## Table of contents

- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Installation](#installation)
- [Examples](#examples)
- [References](#references)

## Introduction

The current approaches formulate mAED problem as a general linear programming (LP) problem [2], which is computationally expensive to solve. This library is proposed for new customized algorithms with three key features: 1) It provides a highly efficient solver to tackle a large and important class of LP problems; 2) It provides a solution for multi-stage decision-making problems with Bayes risk constraints. The package currently supports 2-stage 2-subpopulation AED problem using package [PRIMAL](https://github.com/ShenQianli/primal)[1]. The package will be in future extended to support multi-stage and multi-sub-population AED problem with subspace optimization and row-generation methods[2]. 

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

# User specify D, G, G_c2, loss_table, alpha, beta, J based on data.
est = maed(n_stage=n_stage, n_ppl=n_ppl, prop_ppl=prop_ppl, 
            R_sizes=R_sizes, G_size=G_size,
            D=D, G=G, G_c2=G_c2, alpha=alpha, J=J, beta=beta,
            lambda_0=lambda, loss_tables_0=loss_table,
            lambda_c2=lambda, loss_tables_c2=loss_table,
            solver="primal")

# Visualize the solution path  
plot(est)
```

## References

[1] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.

[2] M. Rosenblum, E. X. Fang, and H. Liu, Optimal, two stage, adaptive enrichment designs for randomized trials, using sparse linear programming, Under minor revision at Journal of Royal Statistical Society, Series B, 2017.


