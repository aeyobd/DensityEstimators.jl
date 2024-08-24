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


## Math-Heavy Background
For those who love ugly equations and generalization.


## Notes


