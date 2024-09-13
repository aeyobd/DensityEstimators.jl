

struct BSplineDensity
    _bspline::BSpline
end


struct LogPSplineDensity
    _bspline::BSpline
end


function log_likelihood(spline, data, weights::Nothing)
    logL = 0.
    for x_i in data
        logL += log(spline(x_i))
    end
    return logL
end


function log_likelihood(spline, data, weights::AbstractVector)
    logL = 0.
    for (x_i, w_i) in zip(data, weights)
        logL += w_i * log(spline(x_i))
    end
    return logL
end


function log_prior(coeffs, prior)
    return log(prior(coeffs))
end


struct VariationalPosterior
    mean::AbstractVector
    log_std::AbstractVector
end

function sample(vp::VariationalPosterior)
    return vp.mean .+ exp.(vp.log_std) .* randn(length(vp.mean))
end


function elbo(vp::VariationalPosterior, data, weights, spline, prior, n_samples)
    elbo = 0.
    for _ in 1:n_samples
        coeffs = sample(vp)
        spline.coeffs .= coeffs
        elbo += log_likelihood(spline, data, weights) + log_prior(coeffs, prior)
    end
    elbo /= n_samples
    return elbo
end
