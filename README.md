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

The current approaches formulate mAED problem as a general linear programming (LP) problem [2], which is computationally expensive to solve. This library is proposed for new customized algorithms with three key features: 1) It provides a highly efficient solver to tackle a large and important class of LP problems; 2) It provides a solution for multi-stage decision-making problems with Bayes risk constraints; 3) It provides additional functions such as visualizing the optimal decision maps. We also provide tutorials on the theoretical background and the code. Please see ``tutorial`` folder for tutorials.

## Directory structure

## Installation

## Examples

## Performance

## References

[1] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.

[2] M. Rosenblum, E. X. Fang, and H. Liu, Optimal, two stage, adaptive enrichment designs for randomized trials, using sparse linear programming, Under minor revision at Journal of Royal Statistical Society, Series B, 2017.


