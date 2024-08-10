import Base: @kwdef


"""
An abstract type for histograms. Should represent a density estimate
sampled at points x with values values and errors err.
"""
abstract type AbstractHistogram end


@kwdef mutable struct Histogram{S, T} <: AbstractHistogram
    bins::AbstractVector{S}
    # x
    # binwidths
    values::AbstractVector{T}
    err::Union{AbstractVector{T}, Nothing} = nothing

    normalization::Symbol = :count
    closed::Symbol = :left
end




@kwdef struct RollingHistogram
    x::AbstractVector
    values::AbstractVector

    bandwidth
    normalization::Symbol
end





"""
    histogram(x[, bins]; weights, normalization, limits, closed)

Computes the histogram of a vector x with respect to the bins with optional weights. Returns the bin edges and the histogram.



Parameters
----------
uncertanties
    :poisson, :bootstrap, :none
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


    hist, low, high = simple_hist(x, bins, weights, closed=closed)

    hist = normalize(hist, bin_volumes(bins), normalization)

    if errors == :poisson
        err = poisson_errors(x, bins, hist)
    else
        err = nothing
    end

    h = Histogram(bins=bins, values=hist, err=err,
                   normalization=normalization, closed=closed)

    return h
end


function poisson_errors(data, bins, hist)
    counts = simple_hist(data, bins, nothing)[1]
    return hist ./ sqrt.(counts)
end


"""
    simple_hist(x, bins[, weights]; closed=:left)

Computes the histogram of a vector x with respect to the bins. Returns the histogram, and the counts for x below and above the histogram range
"""
function simple_hist(x, bins, weights; closed=:left)
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
    idxs = bin_indices(x, bins, closed) .+ 1

    hist = zeros(Int, length(bins)+1)

    for i in idxs
        hist[i] += 1
    end
    
    return hist[2:end-1], hist[1], hist[end]
end



function bin_indices(x, bins, closed=:left)
    if closed == :left
        bin_index = bin_index_left
    elseif closed == :right
        bin_index = bin_index_right
    else
        error("closed must be either :left or :right")
    end

    return [bin_index(bins, x[i]) for i in eachindex(x)]
end



"""
not implemented fully
"""
function bootstrap_hist(x, bins, weights; N=1000)
    bins = make_bins(x, limits, bins)


    for _ in 1:N
        idx = rand(1:length(x), length(x))
        h = simple_hist(x[idx], bins, weights=weights[idx])
        hist += h.values
    end
end


"""
    rolling_histogram(x[, bandwidth]; weights, normalization, limits, samples)

Computes the rolling histogram of a vector x with respect to the bandwidth. Returns the bin edges and the histogram.

"""
function rolling_histogram(x::AbstractVector, bandwidth=bandwidth_freedman_diaconis;
        weights=ones(Int64, length(x)), 
        normalization=:pdf, 
        limits=nothing, 
        samples=10000,
        kwargs...
    )

    if bandwidth isa Function
        bandwidth = bandwidth(x)
    end
    limits = calc_limits(x, limits)

    limits = (limits[1] - bandwidth, limits[2] + bandwidth)

    bins = make_bins(x, limits, samples-1; kwargs...)

    hist = zeros(length(bins))

    N = length(x)

    for i in 1:N
        # range of bins >= x[i] - bandwidth and <= x[i] + bandwidth
        idx_l = searchsortedfirst(bins, x[i] - bandwidth)
        idx_h = searchsortedlast(bins, x[i] + bandwidth)

        if idx_l < 1
            println("idx_l < 1")
        elseif idx_h > length(bins)
            println("idx_h > length(bins)")
        end

        # boundary truncation
        idx_l = max(1, idx_l)
        idx_h = min(length(bins), idx_h)
        width = (idx_h - idx_l + 1)

        hist[idx_l:idx_h] .+= weights[i] / width
    end
    hist = normalize(hist, gradient(bins), normalization)

    h =  RollingHistogram(x=bins, values=hist, bandwidth=bandwidth, normalization=normalization)

    return h
end


function bin_volumes(bins::AbstractVector)
    return diff(bins)
end


"""
Normalizes the values of a histogram.
- :pdf normalizes the histogram to a probability density function
- :count (unnormalized)
- :density normalizes the histogram to the bin width
"""
function normalize(hist::AbstractArray, volumnes::AbstractArray, normalization=:pdf)
    if size(volumnes) != size(hist) 
        error("length of bins must be length of hist + 1")
    end

    if normalization == :pdf
        A = sum(hist) .* volumnes
    elseif normalization == :none
        A = 1
    elseif normalization == :density
        A = volumnes
    elseif normalization == :probabilitymass
        A = sum(hist)
    else
        error("Normalization must be either :none, :pdf, :density or :probabilitymass")
    end

    hist_norm = hist ./ A

    return hist_norm
end



"""
bin_index is index of last bin  <= x
as such bin_index is in 0 to length(bins)
and either extrema represent a value outside the bins
"""
function bin_index_left(bins, x)
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
function bin_index_right(bins, x)
    if x == bins[1]
        return 1
    end
    idx = searchsortedfirst(bins, x) - 1 

    return idx
end




function make_bins(x, limits, bins::Nothing=nothing; bandwidth=nothing)
    if bandwidth == nothing
        error("either bins or bandwidth must be specified")
    end

    bins = limits[1]:bandwidth:(limits[2]+bandwidth)

    return bins
end



function make_bins(x, limits, bins::Int)
    bins = LinRange(limits[1], limits[2], bins+1)

    return bins
end


function make_bins(x, limits, bins::AbstractVector)
    if !issorted(bins)
        error("bins must be sorted")
    end

    if length(unique(bins)) != length(bins)
        error("bins must be unique")
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
        error("bins must be either a number or an array")
    end
end



"""
calculates equal number bins over the array x with n values per bin.
"""
function bins_equal_number(x; n=10)
    return [percentile(x, i) for i in LinRange(0, 100, n+1)]
end


