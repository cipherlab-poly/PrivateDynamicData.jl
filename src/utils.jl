import Distributions
import SpecialFunctions

"""    gaussianMechConstant(ϵ::Real,δ::Real)
Compute the proportionality constant κ(ϵ,δ) for the Gaussian mechanism.
"""
function gaussianMechConstant(ϵ::Real, δ::Real)
    K = sqrt(2) * SpecialFunctions.erfinv(1-2δ)
    return (K+sqrt(K^2+2ϵ))/(2ϵ)
end

"""    (k₁, k₂) = gaussianMechConstant2(ϵ::Real,δ::Real)
Compute the proportionality constant κ(ϵ,δ) for the Gaussian mechanism,
both with the formula using the inverse Q-function (k₁), and according to the
approximation of Theorem A.1 in Dwork and Roth's book (k₂, higher value).
"""
function gaussianMechConstant2(ϵ::Real, δ::Real)
    K = sqrt(2) * SpecialFunctions.erfinv(1-2δ)
    return ((K+sqrt(K^2+2ϵ))/(2ϵ), sqrt(2log(1.25/δ))/ϵ)
end

"""    mean_dp(x::Vector,epsilon::Float64,delta::Float64=0.0)
Compute the average of a set of numbers contained in the vector x,
in an (ϵ,δ)-differentially private way. It is assumed that each entry
x_i of the vector is in [0,1], and the adjacency relation looks at
arbitrary variations in one entry x_i, within this interval.
"""
function mean_dp(x::Vector, ϵ::Real, δ::Real=0.0)
    if δ == 0.0
        d = Distributions.Laplace(0, 1/(ϵ*length(x)))
    else
        d = Distributions.Normal(0, gaussianMechConstant(ϵ,δ)/length(x))
    end
    return Distributions.mean(x) + Distributions.rand(d)
end

"""    laplaceMech(x,f,l1sens::Real,ϵ::Real)
Computes a randomized version to f(x) according to the
Laplace mechanism, for ϵ-differential privacy. l1sens is
the l1-sensitivity of f, which must be computed and provided by the user.
f must take values in R^k, for some k, i.e., return an array of k real values.
"""
function laplaceMech(x, f, l1sens::Real, ϵ::Real)
    t = f(x)
    d = Distributions.Laplace(0, l1sens/ϵ)
    return (t .+ Distributions.rand(d, length(t)))
end

"""    gaussianMech(x,f,l2sens::Real,ϵ::Real,δ::Real)
Computes a randomized version to f(x) according to the
Gaussian mechanism, for (ϵ,δ)-differential privacy. l2sens is
the l2-sensitivity of f, which must be computed by the user.
f must take values in R^k, from some k, i.e., return an array
of k real values.
"""
function gaussianMech(x, f, l2sens::Real, ϵ::Real, δ::Real)
    t = f(x)
    d = Distributions.Normal(0, gaussianMechConstant(ϵ, δ) * l2sens)
    return (t .+ Distributions.rand(d, length(t)))
end

"""    truncatedLaplaceMech(x,f,l1sens::Real,ϵ::Real,δ::Real)
Computes a randomized version to f(x) according to the
truncated Laplace mechanism, for (ϵ,δ)-differential privacy.
This scheme adds bounded noise.
l1sens is the l1-sensitivity of f, which must be computed and provided by the
user. f must take values in R^k, for some k, i.e., return an array of k real
values.

Returns: (r, a) where r is a noisy version of f(x), and a defines the
support of the noise distribution (in the interval [f(x)-a, f(x)+a])
"""
function truncatedLaplaceMech(x, f, l1sens::Real, ϵ::Real, δ::Real)
    t = f(x)
    k = length(t)
    λ = l1sens/ϵ
    a = λ * log( 1 + exp(ϵ) * (k*(1-exp(-ϵ/k)))/(2*δ) )
    # alternative, more conservative but indpt of k
    # a = λ * log(1 + ϵ * exp(ϵ) / (2*δ))
    d1 = Distributions.Laplace(0, λ)
    d = Distributions.truncated(d1, -a, a)
    return (t .+ Distributions.rand(d, length(t)), a)
end
