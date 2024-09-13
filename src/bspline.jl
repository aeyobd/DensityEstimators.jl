

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


"""
    gradient(y[, x])

computes the gradient (dy/dx) of a 2D function at the point (x, y).
assumes that x are sorted.
Returns a vector same length of x with endpoints using linear approximation.
Uses the 2nd order central difference method alla numpy.gradient.
"""
function gradient(y::AbstractVector{T}, x::AbstractVector) where T<:Real
	x = x
	y = y
	N = length(x)

	grad = Vector{T}(undef, N)

	grad[1] = (y[2] - y[1]) / (x[2] - x[1])
	grad[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])
	for i in 2:(N-1)
		hs = x[i] - x[i-1]
		hd = x[i+1] - x[i]

		numerator = hs^2 * y[i+1] + (hd^2 - hs^2) * y[i] - hd^2*y[i-1]
		denom = hd*hs*(hd + hs)
		grad[i] = numerator/denom
	end
	return grad
end



"""
    gradient(y)
computes the gradient
"""
function gradient(y::AbstractVector{T}) where T<:Real
    x = collect(1.0:length(y))
    return gradient(y, x)
end

