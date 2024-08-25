



"""
    binned_statistic(x, y, bins; kwargs...)

Compute a binned statistic for a set of data.

# Keyword Arguments
- `statistic::Function`: The statistic to compute in each bin. Default is `mean`.
- `closed::Symbol`: Which side of the bin is closed. Default is `:left`.
- `weights::Vector`: Weights to apply to the data. Default is `nothing`.
- `limits::Tuple`: The limits of the bins. Default is `nothing`.
"""
function binned_statistic(x, y, bins; 
        statistic=mean,
        closed=:left, 
        weights=nothing,
        limits=nothing
    )

    idxs = bin_indices(x, bins, closed)

    N = length(bins) - 1
    result = Vector{eltype(y)}(undef, N)

    for i in 1:N
        filt = idxs .== i
        if any(filt)
            if weights === nothing
                stat = statistic(y[filt])
            else
                stat = statistic(y[filt], weights[filt])
            end

            result[i] = stat
        else
            result[i] = NaN
        end
    end

    return result
end



"""
    binned_statistic_2d(x, y, z, bins; kwargs...)

Compute a binned statistic for a set of data.
"""
function binned_statistic_2d(x, y, z, bins; 
        statistic=mean,
        closed=:left, 
        weights=nothing,
        limits=nothing
    )

    limits = calc_limits_2d(x, y, limits)
    bins = make_bins_2d(x, y, bins, limits)
    idxs = bin_indices_2d(x, y, bins, closed)
    println(idxs)

    Nx, Ny = length(bins[1]) - 1, length(bins[2]) - 1
    result = Matrix{eltype(z)}(undef, Nx, Ny)

    for i in 1:Nx, j in 1:Ny
        filt = (first.(idxs) .== i) .& (last.(idxs) .== j)
        if any(filt)
            if weights === nothing
                stat = statistic(z[filt])
            else
                stat = statistic(z[filt], weights[filt])
            end

            result[i, j] = stat
        else
            result[i, j] = NaN
        end
    end

    return result
end
