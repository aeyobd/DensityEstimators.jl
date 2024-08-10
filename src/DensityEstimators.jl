module DensityEstimators

using DocStringExtensions: TYPEDEF, FIELDS


export Histogram, Histogram2D


include("interface.jl")
include("bandwidth.jl")
include("bayesian_blocks.jl")
include("histogram.jl")
include("histogram2d.jl")
include("kde.jl")
include("kde2d.jl")
include("knuth_hist.jl")
include("limits.jl")
include("statistics.jl")


end
