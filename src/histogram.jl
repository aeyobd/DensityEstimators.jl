import Base: @kwdef


"""
An abstract type for histograms. Should represent a density estimate
sampled at points x with values values and errors err.
"""
abstract type AbstractDensityEstimate end


"""
A Histogram type.

"""
@kwdef mutable struct Histogram{S, T} <: AbstractDensityEstimate
    bins::AbstractVector{S}
    # x
    # binwidths
    values::AbstractVector{T}
    err::Union{AbstractVector{T}, Nothing} = nothing

    normalization::Symbol = :count
    closed::Symbol = :left
end



"""
    histogram(x[, bins]; weights, normalization, limits, closed)

Computes the histogram of a vector x with respect to the bins with optional weights. Returns the bin edges and the histogram.



# Keyword Arguments
- `weights=nothing`: Weights for each value in x. If nothing, all values are weighted equally.
- `normalization=:none`: Normalization of the histogram. Options are :none, :pdf, :density, :probabilitymass
- `limits=nothing`: The limits of the histogram. If nothing, the limits are the minimum and maximum of x.
- `closed=:left`: The side of the bins that is closed. Options are :left and :right
- `errors=:poisson`: The type of errors to compute. Options are :poisson and :none
- `kwargs...`: Additional keyword arguments to pass to the binning function.



# Examples
```julia
x = [1,1,2,3,4,5]
h = histogram(x, 0:6)
>>> Histogram{Int64, Float64}([0, 1, 2, 3, 4, 5, 6], [0.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0], nothing, :count, :left)
```


"""
function histogram(x::AbstractVector, bins=bandwidth_freedman_diaconis; 
        weights=nothing,
        normalization=:none,
        limits=nothing,
        closed=:left,
        errors=:poisson,
        kwargs...
    )

    limits = calc_limits(x, limits)
    bins = make_bins(x, limits, bins; kwargs...)

    counts, counts_low, counts_high = simple_hist(x, bins, weights, closed=closed)
    hist = normalize(counts, bin_volumes(bins), normalization)
    err = calc_hist_errors(x, bins, hist; errors=errors)

    return Histogram(bins=bins, values=hist, err=err,
        normalization=normalization, 
        closed=closed
    )
end





"""
    simple_hist(x, bins[, weights]; closed=:left)

Computes the histogram of a vector x with respect to the bins. Returns the histogram, and the counts for x below and above the histogram range
"""
function simple_hist(x, bins, weights; closed=:left)
    if length(x) != length(weights)
        throw(DimensionMismatch("length of x and weights must be the same. Got $(length(x)) and $(length(weights))"))
    end

    N = length(bins)
    hist = zeros(eltype(weights), N+1)

    idx = bin_indices(x, bins, closed) .+ 1

    for i in eachindex(x)
        hist[idx[i]] += weights[i]
    end

    return hist[2:end-1], hist[1], hist[end]
end



function simple_hist(x, bins, weights::Nothing=nothing; closed=:left)
    N = length(bins)

    # adding one since we shift array to allow overflow to left and right
    idxs = bin_indices(x, bins, closed) .+ 1

    hist = zeros(Int, length(bins)+1)

    for i in idxs
        hist[i] += 1
    end
    
    return hist[2:end-1], hist[1], hist[end]
end



"""
    bin_indices(x, bins; closed=:left)

Returns the indices of the bins that the values in x fall into.
Closed can be either :left or :right, and determines whether the left or right 
side of the bin is closed.
Values less than the first bin are returned as 0, and values greater than the last bin or NaN are assigned index N.
"""
function bin_indices(x, bins, closed=:left) 
    if length(bins) < 2
        throw(ArgumentError("bins must have at least 2 elements"))
    end

    if closed == :left
        bin_index = _bin_index_left
    elseif closed == :right
        bin_index = _bin_index_right
    else
        throw(ArgumentError("closed must be either :left or :right, got $closed"))
    end

    return [bin_index(bins, x[i]) for i in eachindex(x)]
end



"""
bin_index is index of last bin  <= x
as such bin_index is in 0 to length(bins)
and either extrema represent a value outside the bins
"""
function _bin_index_left(bins, x)
    if x == bins[end]
        return length(bins) - 1
    end
    idx = searchsortedlast(bins, x) 

    return idx
end



"""
bin_index is index of last bin  <= x
as such bin_index is in 0 to length(bins)
and either extrema represent a value outside the bins
"""
function _bin_index_right(bins, x)
    if x == bins[1]
        return 1
    end
    idx = searchsortedfirst(bins, x) - 1 

    return idx
end



"""
    make_bins(x, limits, bins; kwargs...)

Creates the bins for the histogram. 

Bins can be specified as:
- a number of bins
- a vector of bin edges
- a function that returns the bandwidth

Or bins can be left unspecified and the bandwidth may be a number or a function instead.
"""
function make_bins(x, limits, bins::Nothing=nothing; bandwidth=nothing)
    if bandwidth == nothing
        throw(ArgumentError("bins or bandwidth must be specified"))
    end

    bins = limits[1]:bandwidth:(limits[2]+bandwidth)

    return bins
end



function make_bins(x, limits, bins::Int)
    return LinRange(limits[1], limits[2], bins+1)
end


function make_bins(x, limits, bins::AbstractVector)
    if !issorted(bins)
        throw(ArgumentError("bins must be sorted"))
    end

    if length(unique(bins)) != length(bins)
        throw(ArgumentError("bins must be unique"))
    end
    return bins
end


function make_bins(x, limits, bins::Function; kwargs...)
    h = bins(x; kwargs...)
    if h isa Real
        return make_bins(x, limits, bandwidth=h)
    elseif h isa AbstractVector
        return h
    else
        throw(ArgumentError("bins must be a number, a vector or a function that returns a number or a vector"))
    end
end


"""
    bin_volumes(bins)

Returns the volumes of the bins. In 1D, this is the bin width.
Used for normalization
"""
function bin_volumes(bins::AbstractVector{<:Real})
    return diff(bins)
end



"""
    normalize(hist, volumes, normalization=:pdf)

Normalizes the values of a histogram given bin volumes. Valid values of
normalization are:
- :none (unnormalized). Applys no normalization. histogram values represent the
    mass of the data in each bin.
- :pdf normalizes the histogram to a probability density function (i.e. the
    integral of the histogram is 1)
- :density normalizes the histogram to be a mass density (i.e. the integral of
    the histogram is the sum of the values)
- :probabilitymass normalizes the histogram to the sum of the values 
    (i.e. sum of values is unity, not a density)
"""
function normalize(hist::AbstractArray{<:Real}, volumes::AbstractArray{<:Real}, normalization=:pdf)
    if size(volumes) != size(hist) 
        throw(DimensionMismatch("length of bins must be length of hist + 1. Got $(size(volumes)) and $(size(hist))"))
    end

    if normalization == :pdf
        A = sum(hist) .* volumes
    elseif normalization == :none
        A = 1
    elseif normalization == :density
        A = volumes
    elseif normalization == :probabilitymass
        A = sum(hist)
    else
        throw(ArgumentError("Normalization must be either :none, :pdf, :density or :probabilitymass. Got $normalization"))
    end

    hist_norm = hist ./ A

    return hist_norm
end


function calc_hist_errors(x, bins, hist; errors=:poisson)
    if errors == :poisson
        err = poisson_errors(x, bins, hist)
    else
        err = nothing
    end

    return err
end

"""
    poisson_errors(x, bins, hist)

Computes the poisson errors for a histogram given the values of x and the histogram.
"""
function poisson_errors(x, bins, hist)
    counts = simple_hist(x, bins, nothing)[1]
    return hist ./ sqrt.(counts)
end
