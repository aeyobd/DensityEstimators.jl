
struct BayesianBlocks
    edges::Vector{Float64}
    counts::Vector{Int}
end


function p0_prior(N::Int; p0=0.05)
    return 4 - log(73.53 * p0 * (N^-0.478))
end



"""
reproduction of astropy code for bayesian blocks

TODO: refactor, split apart, and validate

"""
function bayesian_blocks(t::Vector{Float64}, x=nothing, fitness::Function=events_fitness;
        prior=p0_prior
    )
    t, x = _validate_input(t, x)
    edges = vcat(t[1], midpoints(t), t[end])


    N = length(t)
    best = zeros(N)
    last = zeros(Int64, N)

    for i in 1:N
        kwargs = Dict{Symbol, Any}()
        kwargs[:N_k] = reverse(cumsum(x[i:-1:1]))
        kwargs[:T_k] = edges[i+1] .- edges[1:i]

        fit_vec = fitness(;kwargs...)

        A_R = fit_vec .- prior(N)
        A_R[2:end] .+= best[1:i-1]

        i_max = argmax(A_R)
        best[i] = A_R[i_max]
        last[i] = i_max
    end

    change_points = zeros(Int64, N)

    i_cp = N+1
    ind = N+1
    while i_cp > 1
        i_cp -= 1

        change_points[i_cp] = ind
        if ind == 1
            break
        end

        ind = last[ind-1]
    end

    if i_cp == 1
        change_points[1] = 1
    end

    change_points = change_points[i_cp:end]

    return edges[change_points]
end


function _validate_input(t::Vector{Float64}, x)
    if isnothing(x)
        x = ones(length(t))
    end
    if length(t) != length(x)
        throw(ArgumentError("t and x must have the same length"))
    end


    filt = sortperm(t)
    return t[filt], x[filt]
end

function events_fitness(; N_k, T_k)
    return @. N_k * log(N_k / T_k)
end

function argnames(f::typeof(events_fitness))
    return (:N_k, :T_k)
end
