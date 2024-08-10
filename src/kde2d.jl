import QuadGK: quadgk



@doc raw"""
    kernel2d_gaussian(x, y)

A 2D gaussian kernel function defined as
```math
f(x, y) = \frac{1}{2\pi} \exp\left(-\frac{x^2 + y^2}{2}\right)
```
"""
function kernel2d_gaussian(x::Real, y::Real)
    return 1/(2π) * exp(-0.5 * (x^2 + y^2))
end



@doc raw"""
    kernel2d_epanechnikov(x, y)

A 2D Epanechnikov kernel function defined as
```math
f(x, y) = \begin{cases}
    \frac{3}{4} (1 - x^2 - y^2) & \text{if } |x| \leq 1 \text{ and } |y| \leq 1 \\
    0 & \text{otherwise}
\end{cases}

See also [`kernel_epanechnikov`](@ref).
"""
function kernel2d_epanechnikov(x::Real, y::Real)
    if x^2 + y^2 > 1
        return 0
    else
        return 2/π * (1 - x^2 - y^2)
    end
end


function integrate2d(kernel::typeof(kernel2d_epanechnikov), x1::Real, x2::Real, y1::Real, y2::Real)

    if x1^2 + y1^2 < 1 && x2^2 + y2^2 < 1
        dx = x2 - x1
        dy = y2 - y1
        dx3 = x2^3/3 - x1^3/3
        dy3 = y2^3/3 - y1^3/3
        return 2/π * (dx * dy - dx3 * dy - dy3 * dx)
    elseif x1^2 + y1^2 < 1
        x2 = sqrt(1 - y2^2)
        dx = x2 - x1
        dx3 = x2^3/3 - x1^3/3
        return 2/π * (dx * y2 - dx3 * y2)
    end
end


"""
A 2D KDE object.
"""
@kwdef mutable struct KDE2D
    x::Vector{F}
    y::Vector{F}
    values::Matrix{F}
    kernel::Function
    r_trunc::F = 3
end




"""
    kde2d(x, y, bandwith; kwargs...)

A 2D kernel density estimator. 

Bandwidth may be 

# Arguments
- `bins=100` The bins to sample the KDE on
- `r_trunc::Real=5.0` The truncation range for the KDE
- `limits=nothing` The limits of the KDE. Expects a tuple of tuples
- `weights=ones(length(x))` The weights for each point
- `kernel=kernel_2d_gaussian` The kernel function to use

"""
function kde2d(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, bandwidth::AbstractVector{<:Real}; 
        bins=100, 
        r_trunc::Float64=5.0, 
        limits=nothing,
        weights=ones(length(x)), 
        kernel::Function=kernel_2d_gaussian
    ) 

    
    xlims, ylims = split_limits(limits)
    xlims = calc_limits(x, xlims)
    ylims = calc_limits(y, ylims)

    xgrid = make_bins(x, xlims, bins-1)
    ygrid = make_bins(y, ylims, bins-1)

    weights = weights / sum(weights)

    kde = KDE2D(x=xgrid, y=ygrid, 
                values=zeros(Float64, length(xgrid), length(ygrid)),
            kernel=kernel, r_trunc=r_trunc)

    # 2D KDE computation by looping over data points
    for k in 1:length(x)
        add_point!(kde, x[k], y[k], bandwidth[k], weights[k])
    end
    
    return kde
end


function kde2d(x, y, bandwidth::Function=bandwidth_knn; kwargs...)
    return kde2d(x, y, bandwidth(x, y); kwargs...)
end


function kde2d(x, y, bandwidth::Real; kwargs...)
    return kde2d(x, y, fill(bandwidth, length(x)); kwargs...)
end



function add_point!(kde::KDE2D, x::Real, y::Real, bandwidth::Tuple, weight::Real=1)
    idx_x, idx_y = kde_grid_indices(kde, x, y, bandwidth)

    for i in idx_x[1]:idx_x[2]
        for j in idx_y[1]:idx_y[2]
            zx = (x - kde.x[i]) / bandwidth[1]
            zy = (y - kde.y[j]) / bandwidth[2]

            kde.values[i, j] += weight * kde.kernel(zx, zy) / (bandwidth[1] * bandwidth[2])
        end
    end
end


function add_point!(kde::KDE2D, x::Real, y::Real, bandwidth::Real, weight::Real=1)
    add_point!(kde, x, y, (bandwidth, bandwidth), weight)
end


function kde_grid_indices(kde::KDE2D, xi::Real, yj::Real, bandwidth::Tuple)
    hx, hy = bandwidth

    dx = kde.r_trunc * hx
    dy = kde.r_trunc * hy
    # Find the indices of the grid points that are within the truncation range
    #
    idx_x_min = bin_index_safe(kde.x, xi - dx)
    idx_x_max = bin_index_safe(kde.x, xi + dx)
    idx_y_min = bin_index_safe(kde.y, yj - dy)
    idx_y_max = bin_index_safe(kde.y, yj + dy)

    return (idx_x_min, idx_x_max), (idx_y_min, idx_y_max)
end


