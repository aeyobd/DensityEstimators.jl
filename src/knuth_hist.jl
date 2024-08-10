import SpecialFunctions: loggamma

"""
    knuth_F(x, M; kwargs...)

Given the number of bins (M), compute the Knuth histogram likelihood for the data x.
"""
function knuth_F(x::AbstractArray, M; kwargs...)
    N = length(x)

    # Compute the histogram
    h = histogram(x, M; kwargs...)
    n_K = h.values

    log_p = (
        N*log(M)
        + loggamma(M/2)
        - M * loggamma(1/2)
        - loggamma(N + M/2)
        + sum(
            loggamma.(n_K .+ 1/2)
        )
        # + K
       )

    return log_p
end



