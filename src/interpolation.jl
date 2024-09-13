import BSplineKit as BSK


"""
    BSpline(knots::AbstractVector, coeffs::AbstractVector, degree::Int)

A basis spline. To evaluate at a point `x`, use `spline(x)`.
Presently a thin wrapper around `BSplinesKit.BSpline`.

# Properties
- `degree::Int`: The degree of the spline.
- `knots::Vector{T}`: The knots of the spline.
- `weights::Vector{T}`: The coefficients of the spline.
"""
struct BSpline
    _spline::BSK.Spline
end


function BSpline(knots::AbstractVector, coeffs::AbstractVector, degree::Int)

    Nc = length(coeffs)
    Nc_exp = length(knots) + degree  - 1
    if Nc != Nc_exp
        throw(DimensionMismatch("Number of coefficients must be equal to the number of knots minus the degree of the spline, got $Nc instead of $Nc_exp"))
    end

    basis = BSK.BSplineBasis(BSK.BSplineOrder(degree + 1), copy(knots))
    spline = BSK.Spline(basis, copy(coeffs))
    return BSpline(spline)
end


function Base.getproperty(spline::BSpline, prop::Symbol) 
    if prop == :degree
        return BSK.order(spline._spline) + 1
    elseif prop == :knots
        return BSK.knots(spline._spline)
    elseif prop == :coefs
        return BSK.coefficients(spline._spline)
    else
        return getfield(spline, prop)
    end
end

function (spline::BSpline)(x)
    return spline._spline(x)
end



"""
    integral(spline::BSpline)

Return the integral of the spline as a new spline
"""
function integral(spline::BSpline)
    spline_int = BSK.integral(spline._spline)

    return BSpline(spline_int)
end



"""
    derivative(spline::BSpline[, N::Int=1])

Return the (N-th) derivative of the spline as a new spline
"""
function derivative(spline::BSpline, N::Int=1)
    spline_diff = BSK.diff(spline._spline, BSK.Derivative(N))

    return BSpline(spline_diff)
end


@doc raw"""
    bspline_basis(knot::Vector{Float64}, i::Int, n::Int, x::Float64)

The basis evaluated at x for the i-th node of a B-spline of order n with knots `knot`.

Here, the B-splines are defined by the Cox-de-Boor recursion formula.

```math
B_{i,0}(x) = 1 \quad \text{if } t_i \leq x < t_{i+1} \quad \text{else } 0
```

```math
B_{i, k}(x) = \frac{x - t_i}{t_{i+k} - t_i} B_{i,k-1}(x) + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}} B_{i+1,k-1}(x)
```

"""
function bspline_basis(knot::AbstractVector{<:Real}, i::Int, n::Int, x::Float64)
    # Base case: 0th order spline (piecewise constant)
    if i > length(knot) - 1 - n
        throw(DomainError("i is too large"))
    elseif i <= 0
        throw(DomainError("i must be positive"))
    end


    if n == 0
        if i == length(knot) - 1 - n
            return (knot[i] <= x <= knot[i+1]) ? 1 : 0
        end
        return (knot[i] <= x < knot[i+1]) ? 1 : 0
    end

    # Recursive definition of B-spline
    term1 = 0.0
    if (knot[i+n] - knot[i]) != 0
        term1 = (x - knot[i]) / (knot[i+n] - knot[i]) * bspline_basis(knot, i, n-1, x)
    end

    term2 = 0.0
    if (knot[i+n+1] - knot[i+1]) != 0
        term2 = (knot[i+n+1] - x) / (knot[i+n+1] - knot[i+1]) * bspline_basis(knot, i+1, n-1, x)
    end

    return term1 + term2
end




function eval_spline(knot::AbstractVector, weights::AbstractVector, n::Integer, x)
    i_low = max(1, searchsortedfirst(knot, x) - 1)
    i_high = min(length(knot)-1-n, i_low + 2n - 1)
    i_low = 1
    i_high = length(knot) - 1 - n

    if x == knot[end]
        return weights[end]
    end
    return sum([bspline_basis(knot, i, n, x) * weights[i] for i in i_low:i_high])
end


function pad_knots(knot::AbstractVector, degree::Int)
    n = length(knot)
    p = degree

    # Pad the knot vector with p additional knots at the beginning and end
    knot_padded = zeros(n + 2p)
    knot_padded[p+1:n+p] .= knot

    # Pad the beginning and end with the first and last knot value
    #
    dx1 = knot[2] - knot[1]
    knot_padded[1:p] .= knot[1] .- dx1 .* (p:-1:1)

    dx2 = knot[n] - knot[n-1]
    knot_padded[n+p+1:n+2p] .= knot[n] .+ dx2 .* (1:p)

    return knot_padded
end


"""
Given a set of data points and a set of knots, 
calculates the spline through integration in the primal basis.
Note that this is not necessarily the most accurate representation
as splines are not orthoganal functions (except for degree 0 == histogram).
"""
function PrimalBSpline(x::AbstractVector, knots; degree=3, normalize=true)
    k1 = pad_knots(knots, degree)

    Nc = length(knots)+degree-1

    a_unnorm = [sum(bspline_basis.([k1], i, degree, x)) for i in 1:Nc]
    a_unnorm ./= length(x)

    spline = BSpline(knots, a_unnorm, degree)

    if normalize
        A = integral(spline)(spline.knots[end])
        a = a_unnorm ./ A
        spline = BSpline(knots, a, degree)
    end


    return spline
end
