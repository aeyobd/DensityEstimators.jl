# DensityEstimators

[![Build Status](https://github.com/aeyobd/DensityEstimators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aeyobd/DensityEstimators.jl/actions/workflows/CI.yml?query=branch%3Amain)

A variety of routines for estimating the probability density given
samples.


## Introduction

How much statistics can you do without looking at some kind of density estimate? How many histograms do you look at throughout your work?
Because this is an ubiquitous problem, I have collected together common and advanced techniques for nonparametric density estimation (i.e. pretty much anywhere where a histogram would be appropriate) under a unified interface to allow for easy switching between methods. 
In this package, I include estimators such as histograms, kernel density estimation. In the figure, I hope to add novel methods including Dirichlet diffusion trees and Dirichlet process regression and monte-carlo and bayesian histogram methods.
Additionally, I include "applied histogram methods" such as binned mean and median trends and am working on expanding to N dimensions and maybe novel metric spaces!

Please beware that this is very much in its infancy and I am just a grad student so I am not going to solve the field of nonparametric density estimation and I would greatly appreciate any help of comments if you somehow encounter this little project and would like to help out.


## Interface

Partly inspired by different interfaces, I am going to build the interface using two steps: model definition and then a model "fitting" to data step. As such, all models will have Struct constructors and then implement a `fit` method taking in the data. The `fit` function then returns a `DensityEstimate` or one of his children to store the model result. Additionally, there is a functional interface so you can use `histogram` and friends without needing two lines.


## Examples


```julia 
julia> x = [1,2,3,2,3,4,3]
julia> histogram(x)

julia> kde(x)


julia> model = Histogram(bins=10)

julia> fit(model, x)
```

## Statistical Background

### Histograms



### Rolling Histograms
A rolling histogram



# Mathematical Background

For those who love ugly equations and generalization.

## Bandwidth methods

For most methods below, there is a fundamental problem that the performance of the method depends on the selection of the bandwidth. As such, countless methods of bandwidth or number of bins have been proposed, ranging from heuristics to simulation-based and bayesian methods. 

### Heuristics

### Likelihood & Cross-validation based

### Bayesian blocks

### Knuth's method



## Histograms

### Classic histograms

A histogram 

### Rolling histograms

Trivial change to the histogram is to use a sliding window rather than a fixed partition of the set. This method has the same number of degrees of freedom as the fix-width histogram (set by the bin-width only) but retains slightly more characteristics of the data and avoids any biases due to where the bins align with the data. Note that this method is equivalent to KDE (below) with a uniform kernel. In essence, a classical histogram is a sparsely sampled KDE estimator with a uniform bin width.

### Bayesian formulation & adding priors to the bins





## Kernel density estimate

Another classic 



### Adaptive kernel estimation

Unlike for histograms, another approach to KDE is to adaptively chose the kernel for each datapoint (or sampled point) to adjust to the density or bias of the data.



## Splines

Another entirely different approach to density estimation is through the use of splines. There are several methods along these lines as discussed in each section below. In particular, the B-Spline density estimation has several different algorithms to estimate the coefficients. 

In all methods, the idea is to fit some spline function (represented by a kernel, a knot vector, and a coefficient vector) which represents a generalized function class to some aspect of the density, such as the CDF, the logarithmic density, or the density directly.  One major advantage of this method over Histograms is that splines represent continuous functions, not step function in the density so are better able to estimate derivatives of the function. Additionally, splines are intuitive, fast methods. A challenge is the statistical theory is younger, and the methods can be much more computationally expensive depending on the method used to calculate the spline coefficients. Ultimanty, the most direct method I know of to robustly estimate the spline coefficient errors is a full MCMC simulation (with dimension equal to the number of knots), which is likely untractable for many datasets, but likely is significantly more robust than a histogram density. 

Note that splines perform the same density reduction as a histogram with the same number of bins, and is a nonparametric density estimate.

### Log density splines

Perhaps the most direct methods of spline fitting simply fit the splines to the logarithm of the density. This avoids a



### P-splines

### Classic B-splines 

The definition of a bspline is 
$$
B_{i, 0}(x) = \begin{cases}
c & x_i \leq x < x_{i+1} \\
0 & {\rm otherwise}
\end{cases}
$$
where $\theta$ is the heaviside step function (1 for nonzero numbers, 0 otherwise). and $c$ is some normalization coefficient.
$$
B_{i, d}(x) = \frac{x-x_i}{x_{i + d} - x_i} B_{i, d-1}(x) + \frac{x_{i+d+1} - x}{x_{i+d+1} - x_{i+1}} B_{i+1, d-1}(x)
$$


The idea with this method is to fit (or estimate) the coefficients of a _B-spline_ to the data using any number of the methods below. So, we would like to project the data onto a b-spline basis with a knot vector $x_k$ where
$$
f(x) \approx \sum_k \alpha_k \phi_{w_k} (x - x_k)
$$
where $\phi_w$ is the transformed spline basis $\phi_w = 1/w\; \phi(x/w)$, $\alpha_k$ is the coefficient for each knot vector $x_k$. 

### Duality basis estimation

### Galerkin Estimator

### Fourier estimation

### Maximization methods





## Dirchlect density estimation and other sophisticated models

