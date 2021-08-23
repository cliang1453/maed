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

Multi-stage Adaptive Enrichment Design (mAED) problem is mostly formulated as a general linear programming (LP) problem [2], which is computationally expensive to solve. This library is proposed for implementation of new customized algorithms with several key features: 1) It provides a highly efficient solver to tackle a large and important class of LP problems; 2) It provides a solution for multi-stage decision-making problems with Bayes risk constraints. This package currently supports solving m-stage 2-subpopulation AED problem, with the core optimization engine supported by [PRIMAL](https://github.com/ShenQianli/primal)[1]. This package provides an interface to the core LP solver, where the users could specify the mAED problem settings, e.g., prior distributions, stage settings, subpopulation settings and additional constraints. The interface formulates the LP problem based on the user specified problem settings, and is implemented based on the Eigen3 library[3]. 

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

If the multi-stage simulated results are available, our program can directly compute the optimal objective value based on such results. we encourage the users to formulate their results into a ```prob_reward.txt``` file. Consider an mAED with m stages and J loss functions, then each row should be formatted as follows: 

```
state_1 action_1 ... state_n action_n trail_probability loss_1 ... loss_J
```
Each entry should be seperated by a tab. The users can then call ```maed``` function as follows:

```R
library(maed)

# See maed.R for detailed description of the input arguments.
est = maed(n_stage=3,
           state_dim=c(4,4,4),
           decision_dim=c(3,3,3),
           prob_reward_file="reward_prob.txt", # this is an example file containing simulated mAED with 3 stages, 3 states and 4 actions.
           J=1,
           beta=c(0.1))
```

We further provide a 2 stages, 2 sub-populations framework for users who wish to integrate the trails simulation into this package (See maed.R for details). When developing upon the provided framework, the users can call ```maed``` function as follows:

```R
library(maed) 

# See maed.R for detailed description of the input arguments.
est = maed(n_stage=2, 
           state_dim=c(4,4), 
           decision_dim=c(3,3), 
           J=2, 
           beta=c(0.1,0.2),
           n_ppl=2, prop_ppl=c(0.5,0.5), 
           reward=list(c(-0.3,-0.1,0.1,0.3), c(-0.3,-0.1,0.1,0.3)),
           lambda_0="normal", loss_tables_0=NULL,
           lambda_c2="normal", loss_tables_c2=NULL,
           G=NULL, G_c2=NULL, G_size=500, alpha=NULL)
```


## What's next
- This package will be in future extended to support advanced methods including row generation and subspace optimization[2], and will be released through CRAN.
- This package will be in future support the visualization of the optimal decision map.

## References

[1] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.  
[2] M. Rosenblum, E. X. Fang, and H. Liu, Optimal, two stage, adaptive enrichment designs for randomized trials, using sparse linear programming, Under minor revision at Journal of Royal Statistical Society, Series B, 2017.  
[3] D. Eddelbuettel and R. François, “Rcpp: Seamless R and C++ integration,” Journal of Statistical Software, vol. 40, no. 8, pp. 1–18, 2011.  
[4] Guennebaud, G., Jacob, B. et al.(2010).  Eigen v3.  http://eigen.tuxfamily.org.  
