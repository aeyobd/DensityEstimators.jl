import StatsBase as sb


function mean(x::AbstractArray; dims=:)
    return sb.mean(x, dims=dims)
end


function mean(x::AbstractArray, w::AbstractArray; dims=:)
    return sb.mean(x, sb.weights(w), dims=dims)
end


"""
Standard error of the mean.
"""
function sem(x::AbstractArray; dims=:)
    return sb.sem(x, dims=dims)
end


"""
Standard deviation.
"""
function std(x::AbstractArray; dims=:)
    return sb.std(x, dims=dims)
end

function std(x::AbstractArray, w::AbstractArray)
    return sb.std(x, sb.weights(w))
end

