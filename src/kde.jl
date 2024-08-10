import QuadGK: quadgk

F = Float64



@doc raw"""
    gaussian_kernel(x)

A gaussian kernel function defined as
```math
f(x) = \frac{1}{\sqrt{2\pi}} \exp\left(-\frac{x^2}{2}\right)
```
"""
function gaussian_kernel(x::Real)
    return 1/sqrt(2π) * exp(-x^2/2)
end


@doc raw"""
    epanechnikov_kernel(x)

An Epanechnikov kernel function defined as
```math
f(x) = \frac{3}{4} (1 - x^2) \quad \text{if} \quad |x| < 1
```
"""
function kernel_epanechnikov(x::Real)
    return 3/4 * (1 - x^2) ? abs(x) < 1 : 0
end


function integrate(k::typeof(kernel_epanechnikov), a, b)
    a = max(-1, min(1, a))
    b = max(-1, min(1, b))
    return (b - a) * (3/4 * (b - a) - 1/4 * (b^3 - a^3))
end


function integrate(k::Function, a, b)
    return quadgk(k, a, b)[1]
end


"""
    normalize_kernel(f, r_trunc)

Normalizes the kernel function `f` over the range `[-r_trunc, r_trunc]`
"""
function normalize_kernel(f::Function, r_trunc::Real)
    return integrate(f, -r_trunc, r_trunc)
end


# structures 
#
#
"""
KDE(x, bandwidth; weights=nothing, kernel=gaussian_kernel, r_trunc=3)

# Fields

$(FIELDS)
"""
@kwdef struct KDE
    """sample points for density"""
    x::Vector{F}

    """the density"""
    values::Vector{F}

    """ the bandwidth of the distribution"""
    bandwidth::Union{Vector{F}, Function, F}

    """ The kernel function. Should take one argument (distance / 
    bandwidth) and return the kernel value. 
    Will be normalized by code."""
    kernel::Function
    r_trunc::F = 3
end






function calc_kde(x::AbstractArray, bandwidth::AbstractArray; 
        weights=nothing, 
        kernel=gaussian_kernel, 
        r_trunc=3, 
        limits=nothing,
        n_samples=1000
    )

    low = calc_limits(x .- bandwidth, limits)[1]
    high = calc_limits(x .+ bandwidth, limits)[2]
    limits = (low, high)

    bins = make_bins(x, limits, n_samples)

    if weights == nothing
        weights = ones(length(x))
    end

    weights = weights / sum(weights)

    N = length(x)
    hist = zeros(length(bins))
    kde = KDE(bins, hist, bandwidth, kernel, r_trunc)

    for i in 1:N
        add_point!(kde, x[i], bandwidth[i], weights[i])
    end

    return kde
end


function calc_kde(x, bandwidth::Function=bandwidth_knn;
        weights=nothing,
        kernel=gaussian_kernel,
        r_trunc=3,
        limits=nothing, 
        n_samples=1000,
        η=1,
        kwargs...)

    bandwidth = η * bandwidth(x; kwargs...)
    return calc_kde(x, bandwidth,
                    weights=weights, 
                    kernel=kernel, 
                    r_trunc=r_trunc, 
                    limits=limits, 
                    n_samples=n_samples, 
                   )
end

function calc_kde(x, bandwidth::Real; kwargs...)
    return calc_kde(x, fill(bandwidth, length(x)), kwargs...)
end


function add_point!(kde::KDE, x, bandwidth, weight)
    dx = kde.r_trunc * bandwidth
    idx_l = bin_index_safe(kde.x, x - dx)
    idx_h = bin_index_safe(kde.x, x + dx)

    dens = kde.kernel.((kde.x[idx_l:idx_h] .- x) ./ bandwidth) ./ bandwidth

    kde.values[idx_l:idx_h] .+= weight .* dens
end


# =============================================================================
# utility functions
# =============================================================================
#
"""
    bin_index_safe(bins, x)

Returns the index of the bin that `x` falls into. 
If `x` is outside the range of `bins`, it returns the closest bin.
"""
function bin_index_safe(bins::Array, x::Real)
    idx = bin_index_left(bins, x)
    return max(1, min(idx, length(bins)))
end


