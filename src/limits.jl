

"""
    calc_limits(x[, limits])

Calculates the limits of x. If limits is a tuple, it should be a 2-tuple of
lower and upper limits. If either limit is nothing, then the maximum/minimum of
x is used instead of that limit. Raises an error if the lower limit is greater
than the upper limit.

Parameters
----------

filt_nan : Bool
    If true, filters out NaN values. Applies to limits as well
filt_inf : Bool
    If true, filters out infinite values. Applies to limits as well


"""
function calc_limits(x, limits::Tuple=(nothing, nothing);
        filt_nan::Bool=true, filt_inf::Bool=false)

    lower, upper = limits

    included = _create_filter(filt_nan, filt_inf)

    x = filter(included, x)


    if !included(lower)
        lower = minimum(x)
    end

    if !included(upper)
        upper = maximum(x)
    end

    if lower > upper
        error("Lower limit is greater than upper limit ($lower > $upper)")
    end

    return lower, upper
end


function calc_limits(x, limits::Nothing; kwargs...)
    return calc_limits(x, (nothing, nothing); kwargs...)
end


function calc_limits_2d(x, y, limits; kwargs...)
    xlimits, ylimits = split_limits(limits)
    xlimits = calc_limits(x, xlimits; kwargs...)
    ylimits = calc_limits(y, ylimits; kwargs...)
    return xlimits, ylimits
end




# 2D limits

"""
    split_limits(limits)

Returns a 2-2-tuple or 4-tuple xy limits as a tuple of 2-tuples.
"""
function split_limits(limits::Tuple{Any, Any, Any, Any})
    return (limits[1:2], limits[3:4])
end


function split_limits(limits::Tuple{Any, Any})
    return limits
end


function split_limits(limits::Nothing)
    return nothing, nothing
end




"""
creates a filter for calc_limits.
"""
function _create_filter(filt_nan::Bool, filt_inf::Bool)
    filters = Any[
        x->!isnothing(x)
    ]

    if filt_nan
        push!(filters, x->!isnan(x))
    end
    if filt_inf
        push!(filters, x->!isinf(x))
    end

    function f_filt(x)
        for f in filters
            if !f(x)
                return false
            end
        end
        return true
    end

    return f_filt
end

