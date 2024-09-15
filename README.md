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

All models considered are framed in the following problem: We have some unknown density function $f$ which we would like to fit some general model to. We ignore parameteric models (where the form of the estimator is assumed to be some kind of limited distribution.) We are then given some sample $X_i$ of $N$ points drawn from the unknown probability distribution function. We then consider a variety of models below which make some assumption about the parameter space of the underlying distribution. 

The goal is to evaluate the performance of models relative to the following criteria

- Computational efficiency. How fast is the model to run and how does this scale with sample size?
- Uncertainty estimation: Is the model able to appropriatly quantify it's uncertainty in the inferred underlying distribution
- Accuracy & convergence: how accurate is the model with respect to a variety of underlying distributions.
- Derivatives: (potentially important) Are the derivatives of the model reasonable as well?
- Simplicity: how sophisticated or straight forward is the model
- Generalizability: does the model generalize to complex or n-dimensions? Can the model appropriately handle weighted samples as well?

I attempt to approach each model with a bayesian framework and compare the more unusual models to more familiar ones (e.g. to KDE and histograms).

While this is very much a work in progress, I how that this summary will be useful to someone one day.



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



### Introduction to B-Splines and their properties

According to @deboor1978, the definition of the $i$th B-spline (basis spline) of order $k$ with knots $\{t_i\}$ is given by
$$
B_{i, k, t}(x) = (t_{i+k} - t_{i})\;[t_i,\cdots,t_{i+k}](t - x)_{+}^{k-1}
$$
Where the brackets denote the divided difference on the function $t\rightarrow (t-x)^{k-1}_+$, which is the leading coefficient of the polynomial of degree $d$ (order of the set minus 1) which agrees with the function at each of the points $\{t_i, \cdots, t_{i+k}\}$. 

The divided difference can be defined recursively as
$$
[t_i,\ldots,t_{i+k}]g(t) = \begin{cases}
\frac{[t_{i},\cdots,t_{r-1},t_{r+1}, \cdots_{i+k}]g(t)
- [t_i,\cdots, t_{s-1},t_{s+1},\cdots,t_{i+k}]g(t)
}{t_{s} - t_r} 
& t_s \ne t_r, t_s, t_r \in \{t_i\ldots t_{i+k}\} \\
g^{(k)}(t_i)/k! & t_i=\ldots = t_{i+k}, g\in C^k
\end{cases}
$$
and $[t_i]g(t) = g(t_i)$ is the case for the zeroth-degree (order 1) polynomial and all $t$ are unique and $t$ is ordered.  This definition works as if $p_n(t)$ is the k-th degree polynomial agreeing at al t less than n, then $p_{n}(t) = p_{n-1}(t) + [t_i,\cdots,t_{i+n}]g(t)\;\Pi_{i<n+1} (t - t_i)$

We note the following properties of the basis:

- $B_{i}$ is zero if $x \notin (t_i, t_{i+d+1})$. 
- $B_i$ is $C^n$ continuous at $x$ where $n = k - \nu$ where $\nu$ is the number of $t_i = x$. 
- As a result, $B_i$ is represented as a degree $d$ piecewise polynomial with breakpoints at each $t_i$.
- $B_i$ goes to zero at either end of the domain
- $\sum B_i(t_i) = 1$ 
- Each $B_i$ are linearly independent, such that $B_i$ forms the basis for the piecewise polynomial space with $d + 1 - \nu$ constraints at each breakpoints $\xi_i$

As such, any spline can be represented as a sum of basis splines, so for any $f \in \mathbb{P}_{k, \xi, \nu}$, 
$$
f(x) = \sum_i \alpha_i B_i(x)
$$
 where $\alpha_i \in \mathbb{R}$.



To compute the values of B-splines more efficiently, we can use the recursive definition:
$$
B_{i, 0}(x) = \begin{cases}
c & x_i \leq x < x_{i+1} \\
0 & {\rm otherwise}
\end{cases}
$$
where $\theta$ is the heaviside step function (1 for nonzero numbers, 0 otherwise). and $c$ is some normalization coefficient.
$$
B_{i, d}(x) = \frac{x-x_i}{x_{i + k-1} - x_i} B_{i, d-1}(x) + \frac{x_{i+k} - x}{x_{i+k} - x_{i+1}} B_{i+1, d-1}(x)
$$


Another major advantage of B-splines is that their integrals and derivatives have closed-form representations as B-splines as well.
$$
\frac{d}{dx} \sum_i \alpha_i B_{i, d}= \sum_i \alpha_i^-\;B_{i, d-1} \\
\alpha_i^- = (k-1)\,\frac{\alpha_i - \alpha_{i-1}}{t_{i+k+1} - t_i}
$$

which follows from the recursive definition of the divided difference and the derivative of $(t-x)^d_+ = -d (t-x)^{d-1}_+$. Since the derivative of a spline is a spline, we can also solve for the antiderivative of a spline:
$$
\int_{t_1}^x \sum_i \alpha_i B_{i, k}(x) dx = \sum_i^{s-1} \alpha^+_i B_{i, k+1}(x)\\
\alpha^+_i = \sum_{j=1}^i \alpha_j \frac{t_{j+k} - t_j}{k}
$$

provided $x \leq t_s$. If $\alpha_i = \delta_{i, l}$, then $\alpha_i^+ = \frac{t_{l+k} - t_l}{k}$ and 
$$
\int B_{l, d}(x) \,dx = \alpha_l^+  \sum_s B_{i, k}(t_s) = \alpha_l^+ = \frac{t_{l + k} - t_l}{k}
$$
So that the total area of a spline is given by
$$
\int_{-\infty}^{\infty} \sum_i \alpha_i B_{i, d}(x)\,dx = \frac{1}{k}\sum_i \alpha_i (t_{i+k} - t_i)
$$
As such, we can create a normalized spline basis
$$
M_{i, d}(x) = \frac{k}{t_{i+k} - t_i} B_{i, d}(x)
$$
such that the area under each $M_i(x)$ integrates to unity, so that a function $\hat{f}$ 
$$
\hat{f} = \sum \beta_i M_{i, d}(x)
$$


defines a pdf if
$$
\beta_i \geq 0, \quad\quad\sum \beta_i = 1
$$
where 
$$
\beta_i = \frac{t_{i+k} - t_i}{k} \alpha_i
$$
We can inforce the summation constraint by only fitting all but one $\alpha_i$ and solving for the last one by subtraction: 
$$
\beta_n = 1 - \sum_i^{n-1} \beta_i
$$


### Closed form B-spline estimators



### Galerkin Estimator

### Fourier estimation

### Maximization methods

To ensure that the spline is normalized, we can maximize over every coefficient except for the last one, which is given by the linear combination of the other coefficients:
$$
\alpha_n = \frac{k}{t_{n+k} - t_n} \left(1 - \sum_i^{n-1} \frac{t_{i+k} - t_i}{k} \alpha_i\right) = k - \sum_i^{n-1} \frac{t_{i+k} - t_i}{t_{n+k} - t_n} \alpha_i
$$
so that defining 
$$
c_i \equiv\frac{t_{i+k} - t_i}{t_{n+k} - t_n}
$$
lets us define
$$
\alpha_n = k - \sum_i^{n-1} c_i\, \alpha_i
$$

$$
\frac{\partial \alpha_n}{\partial \alpha_i} = -c_i
$$
and all second derivatives are zero.

As the log likelihood is  $ {\cal L} = \sum_i \log f(X_i)$, we have
$$
{\cal L} = \sum_{x\in X} \log\left(\sum_i \alpha_i B_{i}(x)\right)
$$
so that we know the derivatives
$$
\frac{\partial {\cal L}}{\partial \alpha_i} = \sum_x \frac{1}{f(x)} \left(B_i(x) - c_i B_n(x)\right)
$$

$$
\frac{\partial^2 {\cal L}}{\partial \alpha_i \partial \alpha_j} = \sum_x \frac{-1}{f(x)^2} \left(B_i(x) - c_i\, B_n(x)\right) \left(B_j(x) - c_i\,B_n(x)\right)
$$

where the second term arrises from the normalization constraint on the coefficients. This allows for efficient optimization of the coefficients, and we can furthermore add the constraint that all $\alpha_i \geq 0$ (which includes the last coefficient above). 

### Log density splines

Perhaps the most direct methods of spline fitting simply fit the splines to the logarithm of the density. This avoids a



## Orthoganal series estimation



## Wavelet estimation



## Dirchlect density estimation 







## Bandwidth methods

For almost any method above, a fundamental parameter that needs to be estimated is the bandwidth/number of bins/where to place breakpoints.

As such, countless methods of bandwidth or number of bins have been proposed, ranging from heuristics to simulation-based and bayesian methods. 

### Heuristic methods

### Likelihood & Cross-validation based

### Bayesian blocks

### Knuth's method
