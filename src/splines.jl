import BSplineKit as BSK
# BSK does not behave on boundaries ;-;


"""
    BSpline(knots::AbstractVector, coeffs::AbstractVector, order::Int)

A basis spline created with knots `knots`, coefficients `coeffs` and order `order`.
Presently a thin wrapper around `BSplinesKit.BSpline`.
"""
struct BSpline
    knots::AbstractVector
    coefficients::AbstractVector
    order::Int

    function BSpline(knots::AbstractVector, coeffs::AbstractVector, order::Int)
        Nc = length(coeffs)
        Nc_exp = length(knots) - order

        if Nc != Nc_exp
            throw(DimensionMismatch("Number of coefficients must be the number of knots minus the order. Got $Nc instead of $Nc_exp"))
        end

        if !issorted(knots)
            throw(ArgumentError("Knots must be sorted"))
        end

        return new(knots, coeffs, order)
    end
end



"""
returns the order of the spline
"""
function order_of(spline::BSpline)
    return spline.order
end

"""
returns the knots of the spline
"""
function knots_of(spline::BSpline)
    return spline.knots
end

"""
returns the coefficients of the spline
"""
function coefficients_of(spline::BSpline)
    return spline.coefficients
end


"""
evaluates the spline at x
"""
function (spline::BSpline)(x)
    return eval_spline(knots_of(spline), coefficients_of(spline), order_of(spline), x)
end




"""
    integral(spline::BSpline)

Return the integral of the spline as a new spline
"""
function integral(spline::BSpline)
    # spline_int = BSK.integral(spline._spline)

    knots = knots_of(spline)
    order = order_of(spline)
    coeffs = get_int_coefs(spline)

    return BSpline([knots[1]; knots], coeffs, order+1)
end


function get_int_coefs(spline)
    t = spline.knots
    a = spline.coefs
    d =spline.degree
    
    S = length(a)
    
    alpha=0
    coefs = zeros(S+1)
    
    for i in 1:S
        alpha += a[i] * (t[i + d]  - t[i]) / d
        
        coefs[i+1] = alpha
    end
    
    return coefs
end


function area_of(spline::BSpline)
    t = knots_of(spline)
    α = coefficients_of(spline)
    k = order_of(spline)

    n = length(α)
    return sum((t[i+k] - t[i]) * α[i] for i in 1:n) / k
end



"""
    derivative(spline::BSpline[, N::Int=1])

Return the (N-th) derivative of the spline as a new spline
"""
function derivative(spline::BSpline, N=1)
    t = knots_of(spline)
    k = order_of(spline)

    alpha = coefficients_of(spline)

    alpha_1 = (k-1) * [ (alpha[i] - alpha[i-1]) / (t[i+k-1] - t[i]) for i in 2:length(alpha)]

    spline = BSpline(t[2:end-1], alpha_1, k-1)

    if N == 1
        return spline
    else
        return derivative(spline, N-1)
    end
end


"""
    normalized_spline(knots::AbstractVector, coeffs::AbstractVector, order::Int)

Construct a normalized spline from the given knots and coefficients.
Takes all the spline coefficients except the last one, which is then calculated
to be such that the integrated area under the spline is one.
"""
function normalized_spline(knots::AbstractVector, coeffs::AbstractVector, order::Integer)
    β = [coeffs[i] * (knots[order + i] - knots[i]) for i in 1:length(coeffs)] / order
    β_end = 1 - sum(β)
    α_end = order / (knots[end] - knots[end - order]) * β_end

    spline =  BSpline(knots, [coeffs; α_end], order)

    return spline
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
function bspline_basis(knot::AbstractVector{<:Real}, i::Int, order::Int, x::Float64)
    # Base case: 0th order spline (piecewise constant)
    Nk = length(knot)
    if i > Nk - order
        throw(DomainError("i is larger than the knots minus the order, got $i instead of $(Nk - order)"))
    elseif i <= 0
        throw(DomainError("i must be positive, got $i"))
    end

    if x < knot[i] || x >= knot[i+order]
        return 0.
    end

    if order == 1
        return 1.
    end

    # Recursive definition of B-spline
    term1 = 0.0
    if (knot[i+order-1] - knot[i]) != 0
        term1 = (x - knot[i]) / (knot[i+order-1] - knot[i]) * bspline_basis(knot, i, order-1, x)
    end

    term2 = 0.0
    if (knot[i+order] - knot[i+1]) != 0
        term2 = (knot[i+order] - x) / (knot[i+order] - knot[i+1]) * bspline_basis(knot, i+1, order-1, x)
    end

    return term1 + term2
end




function eval_spline(knot::AbstractVector, weights::AbstractVector, k::Integer, x::Real)
    i_low = max(1, searchsortedfirst(knot, x) - k)
    i_high = min(length(knot)-k, i_low + k)

    if x == knot[end]
        return weights[end]
    end
    return sum([bspline_basis(knot, i, k, x) * weights[i] for i in i_low:i_high])
end

function eval_spline(knot, weights, k, x::AbstractVector)
    return [eval_spline(knot, weights, k, x_i) for x_i in x]
end



"""
Given a knot vector and an order, pads the knot vector with p additional knots at the beginning and end.
"""
function pad_knots(knot::AbstractVector, order::Int)
    n = length(knot)
    p = order - 2

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
function PrimalBSpline(x::AbstractVector, knots; order, normalize=true)
    Nc = length(knots) - order 

    a_unnorm = [sum(bspline_basis.([knots], i, order, x)) for i in 1:Nc]
    a_unnorm ./= length(x)


    if normalize
        A = area_of(BSpline(knots, a_unnorm, order))
        a = a_unnorm ./ A
    else
        a = a_unnorm
    end

    spline = BSpline(knots, a, order)
    return spline
end
