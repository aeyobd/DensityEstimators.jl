import QuadGK: quadgk

F = Float64



@doc raw"""
    gaussian_kernel(x)

A gaussian kernel function defined as
```math
f(x) = \frac{1}{\sqrt{2\pi}} \exp\left(-\frac{x^2}{2}\right)
```
"""
function kernel_gaussian(x::Real)
    return 1/sqrt(2π) * exp(-x^2/2)
end


function support_of(k::typeof(kernel_gaussian))
    return (-Inf, Inf)
end


@doc raw"""
    kernel_epanechnikov(x)

An Epanechnikov kernel function defined as
```math
f(x) = \frac{3}{4} (1 - x^2) \quad \text{if} \quad |x| < 1
```

This kernel is 2nd degree and compact.
"""
function kernel_epanechnikov(x::Real)
    if abs(x) <= 1
        return 3/4 * (1 - x^2)
    else
        return zero(x)
    end
end

function support_of(k::typeof(kernel_epanechnikov))
    return (-1, 1)
end


function L2_kernel_quadratic(x::Real)
    if abs(x) <= 1
        return 15/4 * (3x^2 - 1)
    else
        return zero(x)
    end
end


@doc raw"""
    L2_kernel_quartic(x)

A 4th degree kernel function to estimate the second derivative of 
the density (2nd order bias). Defined as
```math
f(x) = \frac{105}{16} (6x^2 - 5x^4 - 1) \quad \text{if} \quad |x| < 1
```
"""
function L2_kernel_quartic(x::Real)
    if abs(x) <= 1
        return 105/16 * (6x^2 - 5x^4 - 1)
    else
        return zero(x)
    end
end



function integrate(k::typeof(kernel_epanechnikov), a, b)
    a = max(-1, min(1, a))
    b = max(-1, min(1, b))
    return (3/4 * (b - a) - 1/4 * (b^3 - a^3))
end


function d_kernel_epanechnikov(x::Real)
    if abs(x) <= 1
        return -3/2 * x
    else
        return 0
    end
end


function d2_kernel_epanechnikov(x::Real)
    if abs(x) <= 1
        return -3/2
    else
        return 0
    end
end

function differentiate(k::typeof(kernel_epanechnikov))
    return d_kernel_epanechnikov
end


function differentiate(k::typeof(d_kernel_epanechnikov), i::Int)
    if i == 1
        return d_kernel_epanechnikov
    elseif i == 2
        return d2_kernel_epanechnikov
    elseif i > 2
        return x->0
    else
        raise(ArgumentError("i must be a non-negative integer"))
    end
end

function differentiate(k::typeof(d_kernel_epanechnikov), i::Int, x::Real)
    if i == 1
        return d2_kernel_epanechnikov(x)
    elseif i > 1
        return x -> 0
    else
        raise(ArgumentError("i must be a non-negative integer"))
    end
end

function moment_of(k::typeof(kernel_epanechnikov), i::Int)
    if i == 0
        return 1
    elseif i == 1
        return 0
    elseif i == 2
        return 1/5
    else
        return integrate(x -> k(x) * x^i, -1, 1)
    end
end

function order_of(k::typeof(kernel_epanechnikov))
    return 2
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

"""
    RBCKernel(kernel, correction, ρ)

A robust bias corrected kernel. The kernel is defined as
```math
f(x) = K(x) + \rho C(x)
```
"""
@kwdef struct RBCKernel
    kernel::Function
    correction::Function
    ρ::Real = 1
    order::Int = order_of(kernel)
end

function (k::RBCKernel)(x::Real)
    return k.kernel(x) - k.ρ^k.order * k.correction(x * k.ρ) * moment_of(k.kernel, k.order)
end
# structures 
#
#
#


"""
    KDE
"""
@kwdef struct KDE
    bandwidth
    kernel::Function
    r_trunc::F = 3
end


"""
    KDEResult

# Fields

$(FIELDS)
"""
@kwdef struct KDEResult
    """sample points for density"""
    x::Vector{F}

    """the density"""
    values::Vector{F}

    """confidence intervals"""
    ci_low::Vector{F}
    ci_high::Vector{F}

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
        kernel=kernel_gaussian, 
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
    kde = KDEResult(bins, hist, bandwidth, kernel, r_trunc)

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
    return calc_kde(x, fill(bandwidth, length(x)); kwargs...)
end


function add_point!(kde::KDEResult, x, bandwidth, weight)
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
    idx = _bin_index_left(bins, x)
    return max(1, min(idx, length(bins)))
end


