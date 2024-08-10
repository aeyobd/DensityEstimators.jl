@kwdef struct Histogram2D{A, B, T}  <: AbstractHistogram
    xbins::AbstractVector{A}
    ybins::AbstractVector{B}
    values::AbstractMatrix{T}
    outside::Real = 0
end


"""
    histogram2d(x, y, bins; weights, limits) -> Histogram2D

Compute a 2D histogram of the data `x` and `y` with the number of bins in each direction given by `bins`.
If bins is

Parameters
----------
x : AbstractVector
    The x data
y : AbstractVector
    The y data
bins :
    The bin edges in each direction
    - AbstractVector: the bin edges in each direction
    - Tuple{AbstractVector, AbstractVector}: the bin edges in each direction
    - Int: the number of bins in each direction
    - Tuple{Int, Int}: the number of bins in each direction
    - Function: the bin edges in each direction. The function f called 
    with f(x, y; kwargs...).
weights : AbstractVector (optional)
    The weights of each data point
limits : Tuple{Tuple{Real, Real}, Tuple{Real, Real}}
    If bins is an Int, then the limits of the data,
    otherwise ignored
normalization: Symbol
    The normazation should be one of
    - :none: no normalization (counts / sum of weights per bin)
    - :pdf: normalize to a probability density function (counts / sum of weights per bin / bin volume)
    - :probabilitymass: normalize to a probability mass function (counts / sum of weights)
    - :density: normalize to a density (counts / sum of weights / bin volume)

"""
function histogram2d(x::AbstractVector, y::AbstractVector, bins=20; 
        weights=nothing, limits=nothing, normalization=:none,
        kwargs...)::Histogram2D

    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length. Got $(length(x)) and $(length(y))."))
    end


    limits = calc_limits_2d(x, y, limits; kwargs...)
    xbins, ybins = make_bins_2d(x, y, bins, limits)

    H, outside = hist2d_simple(x, y, (xbins, ybins), weights)
    H = normalize(H, bin_volumes(xbins, ybins), normalization)

    return Histogram2D(
        xbins=xbins,
        ybins=ybins,
        values=H,
        outside=outside ./ sum(H)
    )
end


function make_bins_2d(x, y, bins::Nothing, limits)
    return make_bins(x, limits[1], bins), make_bins(y, limits[2], bins)
end

function make_bins_2d(x, y, bins::AbstractVector, limits)
    return bins, bins
end


function make_bins_2d(x, y, bins::Tuple{Any, Any}, limits)
    xbins, ybins = bins
    return make_bins(x, limits[1], xbins), make_bins(y, limits[2], ybins)
end

function make_bins_2d(x, y, bins::Int, limits)
    return make_bins(x, limits[1], bins), make_bins(y, limits[2], bins)
end



"""
    hist2d_simple(x, y, bins[, weights]; closed=:left)

Compute a 2D weighted histogram of the data `x` and `y` with the number of bins in each direction given by `bins`.
"""
function hist2d_simple(x, y, bins, weights; closed=:left)
    idxs = bin_indices_2d(x, y, bins, closed)

    outside = zero(eltype(weights))
    Nx, Ny = length(bins[1]) - 1, length(bins[2]) - 1
    H = zeros(eltype(weights), Nx, Ny)

    for (w, (i, j)) in zip(weights, idxs)
        if 1 <= i <= Nx && 1 <= j <= Ny
            H[i, j] += w
        else
            outside += w
        end
    end

    return H, outside
end


function hist2d_simple(x, y, bins, weights::Nothing=nothing; closed=:left)
    idxs = bin_indices_2d(x, y, bins, closed)

    outside = 0
    Nx, Ny = length(bins[1]) - 1, length(bins[2]) - 1

    H = zeros(Int, Nx, Ny)
    for (i, j) in idxs
        if 1 <= i <= Nx && 1 <= j <= Ny
            H[i, j] += 1
        else
            outside += 1
        end
    end

    return H, outside
end


function bin_indices_2d(x, y, bins, closed=:left)
    xbins, ybins = bins
    if closed == :left
        bin_index = bin_index_left
    elseif closed == :right
        bin_index = bin_index_right
    else
        throw(ArgumentError("Unknown closed method: $closed"))
    end

    return (
        (bin_index(xbins, xi), bin_index(ybins, yi)) 
            for (xi, yi) in zip(x, y))
end


"""
The centres of each bin in the histogram as an array of x and y
"""
function bin_centres(h::Histogram2D)
    xcentres = midpoints(h.xbins)
    ycentres = midpoints(h.ybins)'
    return xcentres .+ 0'ycentres, ycentres .+ 0*xcentres
end



function bin_volumes(xbins, ybins)
    dx = diff(xbins) 
    dy = diff(ybins)
    return dx .* dy' # x -> row num = column vector, y -> column num = row vector
end



