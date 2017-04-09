using Distributions

"""    gaussianMechConstant(ϵ::Real,δ::Real)
Compute the proportionality constant κ(ϵ,δ) for the Gaussian mechanism.
"""
function gaussianMechConstant(ϵ::Real, δ::Real)
	K = sqrt(2) * erfinv(1-2δ)
	return (K+sqrt(K^2+2ϵ))/(2ϵ)
end

"""    (k₁, k₂) = gaussianMechConstant2(ϵ::Real,δ::Real)
Compute the proportionality constant κ(ϵ,δ) for the Gaussian mechanism,
both with the formula using the inverse Q-function (k₁), and according to the
approximation of Theorem A.1 in Dwork and Roth's book (k₂, higher value).
"""
function gaussianMechConstant2(ϵ::Real, δ::Real)
	K = sqrt(2) * erfinv(1-2δ)
	return ((K+sqrt(K^2+2ϵ))/(2ϵ), sqrt(2log(1.25/δ))/ϵ)
end

"""    mean_dp(x::Vector,epsilon::Float64,delta::Float64=0.0)
Compute the average of a set of numbers contained in the vector x,
in an (ϵ,δ)-differentially private way. It is assumed that each entry
x_i of the vector is in [0,1], and the adjacency relation looks at
arbitrary variations in one entry x_i, within this interval.
"""
function mean_dp(x::Vector, ϵ::Real, δ::Real=0.0)
	if δ==0
		d = Laplace(0,1/(ϵ*length(x)))
	else
		d = Normal(0,1)
	end
	return mean(x)+rand(d,1)
end

"""    laplaceMech(x,f,l1Sens::Real,ϵ::Real)
Computes a randomized version to f(x) according to the
Laplace mechanism, for ϵ-differential privacy. l1Sens is
the l1-sensitivity of f, which must be computed by the user.
f must take values in R^k, from some k, i.e., return an array
of k real values.
"""
function laplaceMech(x, f, l1sens::Real, ϵ::Real)
	t = f(x)
	d = Laplace(0,l1sens/ϵ)
	return (t + rand(d,length(t)))
end

"""    gaussianMech(x,f,l2Sens::Real,ϵ::Real,δ::Real)
Computes a randomized version to f(x) according to the
Gaussian mechanism, for (ϵ,δ)-differential privacy. l2Sens is
the l2-sensitivity of f, which must be computed by the user.
f must take values in R^k, from some k, i.e., return an array
of k real values.
"""
function gaussianMech(x, f, l2Sens::Real, ϵ::Real, δ::Real)
	t = f(x)
	d = Normal(0, gaussianMechConstant(ϵ,δ)*l2sens)
	return (t + rand(d,length(t)))
end
