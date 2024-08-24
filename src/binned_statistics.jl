



"""
    binned_statistic(x, y, bins, statistic)

Compute a binned statistic for a set of data.
"""
function binned_statistic(x, y, bins, statistic; 
        closed=:left, 
        weights=nothing,
        limits=nothing
    )

    idxs = bin_indices(x, bins, closed)

    for i in eachindex(y)
        idx = idxs[y]

    end
end



function binned_statistic(x, bins, statistic)
end
