
@kwdef struct RollingHistogram
    x::AbstractVector
    values::AbstractVector

    bandwidth
    normalization::Symbol
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


