# code to interface with varius external packages

import NearestNeighbors as nn
import StatsBase: quantile, percentile



"""
    kth_nns(x, k=5, metric=:euclidean)

Compute the k-th nearest neighbor of each point in `x` using the Euclidean distance.

If x is a NxD matrix, then 
"""
function knn(x::AbstractArray; k=5, metric=:euclidean)
    if x isa AbstractVector
        x = x'
    end
    tree = nn.KDTree(x)
    idxs, dists = nn.knn(tree, x, k+1, true)

    idxs = [d[end] for d in idxs]
    dists = [d[end] for d in dists]
    return idxs, dists
end


