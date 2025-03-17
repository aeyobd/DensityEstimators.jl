# Mathematical Background

For those who love ugly equations and generalization.

All models considered are framed in the following problem: We have some unknown density function $f$ which we would like to fit some general model to. We ignore parameteric models (where the form of the estimator is assumed to be some kind of limited distribution.) We are then given some sample $X_i$ of $N$ points drawn from the unknown probability distribution function. We then consider a variety of models below which make some assumption about the parameter space of the underlying distribution. 

The goal is to evaluate the performance of models relative to the following criteria

- **Computational efficiency**. How fast is the model to run and how does this scale with sample size? (Histograms with poisson errors are an excellent baseline)
- **Accuracy & convergence**. how accurate is the model with respect to a variety of underlying distributions.
- **Uncertainty estimation**. Is the model able to appropriately quantify its uncertainty in the inferred underlying distribution? Do confidence intervals have appropriate coverage properties?
- **Derivatives**. (Potentially important.) Are the derivatives of the model reasonable as well? Does uncertainty quantification adequately cover derivatives as well?
- **Simplicity**. How sophisticated or straight forward is the model?
- **Generalizability**. Does the model generalize to complex or n-dimensions? Can the model appropriately handle weighted samples as well?

I attempt to approach each model with a bayesian framework and compare the more unusual models to more familiar ones (e.g. to KDE and histograms).

While this is very much a work in progress, I how that this summary will be useful to someone one day.

## Terminology / notation

- Vectors are bf, ${\bf b}$. 
- $f$ is the true density and $\hat f$ is the density estimator.
-  $x_i$ is the vector of observed sample, $w_i$ and $\sigma_i$ are the weights and gaussian uncertainties of the samples (if relevant).
- $\sum$s are taken over the length of model coefficients $K$
- $k_i$ is the knot vector / bin edges if relevant
- $\lambda$ is a smoothing parameter
- $\kappa_i$ is the ith kernel (based on bandwidth/knots) and $\kappa$ is the unscaled/normalized kernel generator

# Histograms

### Classic histograms

A histogram is, in essence, a piecewise constant density estimate
$$
\hat f(x) = \sum_i b_i\ \kappa_i(x)\\
\kappa_i(x) = \begin{cases}
1 / (k_{i+1} - k_i) & k_i < x < k_{i+1}\\
0 & {\rm otherwise}
\end{cases}
$$
The maximum likelihood estimate solution for $b_i$ is simply the number of observations in that bin over the total number of observations 
$$
b_i = \frac{N_i}{N}
$$
This gives the familiar representation of a histogram as we know it. 

Similarly, under a frequentist perspective, each bin is represented as a poisson process, resulting in an standard uncertaninty of $\sigma b_i = \sqrt{N_i} / N$. In detail, we can construct confidence intervals from this expression, but often assuming the limiting gaussian case is adequate if there are more than about 5 observations in the bin.
$$
b_i \sim {\rm Poisson}(N_i/N)
$$


Under a bayesian framework, a reasonable prior on the histogram coefficients $b_i$ is a Dirichlet prior
$$
{\bf b} \sim {\rm Dir()}
$$


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

Another framework is to use orthoganal series / wavelets as function expansions for the density. In some sense, splines and histograms are related or special cases of this section. 

## Dirchlect density estimation 







# Bandwidth and Knot selection methods

For almost any method above, a fundamental parameter that needs to be estimated is the bandwidth/number of bins/where to place breakpoints.

As such, countless methods of bandwidth or number of bins have been proposed, ranging from heuristics to simulation-based and bayesian methods. 

### Heuristic methods

### Likelihood & Cross-validation based

### Bayesian blocks

### Knuth's method

### Hierarchical bayesian framework

Perhaps the most robust method of including a bandwidth or knot vector in the model is to simply add it as part of a hierarchical bayesian model. For example, for the bandwidth/$K$/bin width, we can add e.g. an inverse gamma prior to this model and 
$$
x \sim {\rm model}(x_i, w) \\
w \sim {\rm InvGamma}(\alpha, \beta)
$$
