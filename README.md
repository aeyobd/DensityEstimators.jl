# DensityEstimators

[![Build Status](https://github.com/aeyobd/DensityEstimators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aeyobd/DensityEstimators.jl/actions/workflows/CI.yml?query=branch%3Amain)

A variety of routines for estimating the probability density given
samples.


## Introduction

A classic problem in statistics is to estimate an underlying probability
density given some sample of the distribution. Perhaps the most straightforward
method is a simple histogram (implemented here as well), however, given the 
nature of the specific problem, other alternatives to histograms may be 
preferable. For example, Kernel Density Estimation (KDE) and variants thereof,
and there are a number of yet more-sophisticated methods (include Dirichlet 
processes, )
Additionally, there are countless heuristics and models to estimate the bin-size 
or bandwidth of a PDF estimate; so here we implement many of these variants as well.

## Interface


## Examples


## Statistical Background

### Histograms



### Rolling Histograms
A rolling histogram


## Math-Heavy Background
For those who love ugly equations and generalization.


## Notes


