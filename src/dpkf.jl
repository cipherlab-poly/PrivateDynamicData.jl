import Pkg

#=
if haskey(Pkg.installed(),"MosekTools")
    using MosekTools
    const dpkf_ok = true
else
    dpkf_ok = false
    println("WARNING: MosekTools.jl not installed, the functions
    for differentially private Kalman filtering cannot be used!")
end
=#

if haskey(Pkg.installed(),"COSMO")
    using COSMO
    const dpkf_ok = true
else
    dpkf_ok = false
    println("WARNING: COSMO.jl not installed, the functions
    for differentially private Kalman filtering cannot be used!")
end

using JuMP
using ControlSystems
using LinearAlgebra


"""    (D, P_val, X_val, Ω_val, M) =
            staticInputBlock_DPKF_ss(Ls, As, Cs, Winvs, Vs, Vinvs, k_priv, ρ,
                                    m=size(As,1), p=size(Cs,1), r=size(Ls,1))

Compute an input matrix D for the two-block differentially private
steady-state Kalman filter mechanism. D takes linear combinations of
individual signals before privacy-preserving Gaussian noise injection.

Note that we do not guarantee observability of the pair (A,D*C). If
observability is lost and needed, call the dfactor function with the
returned matrix M and a truncation threshold lower than 1e-6 in order
to compute a new matrix D with more rows.

Inputs: one matrix per individual

- Ls: size(r,nusers*m), or size (r,m,nusers) if want estimator of sum_i L_i x_i
- As: size (m,m,nusers)
- Cs: size (p,m,nusers)
- Winvs: size (m,m,nusers) -- inverses of process noise cov. matrices
- Vs: size (p,p,nusers) -- measurement noise cov. matrices
- Vinvs: inverses of Vs (computed externally for modularity/efficiency)
- ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.
"""

function staticInputBlock_DPKF_ss(Ls, As, Cs, Winvs, Vs, Vinvs, ρ, k_priv,
                                  nusers=size(As,3), m=size(As,1),
                                  p=size(Cs,1), r=size(Ls,1))

    #modl = Model(Mosek.Optimizer)
    modl = Model(COSMO.Optimizer, "max_iter" => 20000)
    #modl = Model(with_optimizer(Mosek.Optimizer))
    #modl = Model(with_optimizer(SCS.Optimizer))

    Nm = nusers * m
    Np = nusers * p

    A = zeros(Nm, Nm)
    C = zeros(Np, Nm)
    V = zeros(Np, Np)
    Winv = zeros(Nm, Nm)
    Vinv = zeros(Np, Np)

    for i=1:nusers
        A[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = As[:,:,i]
        C[(i-1)*p+1:i*p, (i-1)*m+1:i*m] = Cs[:,:,i]
        V[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vs[:,:,i]
        Winv[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Winvs[:,:,i]
        Vinv[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vinvs[:,:,i]
    end

    if size(Ls,3) > 1
        L = zeros(r, Nm)
        for i=1:nusers
            L[:, (i-1)*m+1:i*m] = Ls[:,:,i]
        end
    else
        L = copy(Ls)
    end

    if r == 1
        @variable(modl, X >= 0)
        @objective(modl, Min, X)
    else
        @variable(modl, X[1:r, 1:r], PSD)
        @objective(modl, Min, tr(X))
    end
    if Nm == 1
        @variable(modl, Ω >= 0)
    else
        @variable(modl, Ω[1:Nm, 1:Nm], PSD)
    end
    if Np == 1
        @variable(modl, P >= 0)
    else
        @variable(modl, P[1:Np, 1:Np], PSD)
    end

    @SDconstraint(modl, hvcat((2,2), X, L, L', Ω) >= zeros(r+Nm,r+Nm))

    @SDconstraint(modl, hvcat((2,2),
        C'*P*C-Ω+Winv, Winv*A, A'*Winv, Ω+A'*Winv*A) >= zeros(2*Nm,2*Nm))

    for i=1:nusers
        e = zeros(p,Np); e[:,(i-1)*p+1:i*p]=Matrix{Float64}(I,p,p)
        @SDconstraint(modl,
          hvcat((2,2), Matrix{Float64}(I,p,p)/(k_priv^2*ρ[i]^2)+Vinvs[:,:,i], e, e', V-V*P*V)
          >= zeros(Np+p,Np+p))
    end

    optimize!(modl)
    status = termination_status(modl)

    println("============================================")
    println("Solution status: ", status)
    println("Objective value: ", objective_value(modl))
    println("============================================")

    P_val = value.(P)

    X_val = value.(X)
    Ω_val = value.(Ω)

    #= # Debug
    Σ_val = inv(Ω_val)
    println("Norm of the difference between Σ^{-1} and (A Σ A' + W)^{-1}+C' P C")
    E1 = C'*P_val*C+inv(A*Σ_val*A'+inv(Winv))-inv(Σ_val)
    println(norm(E1))
    =#

    # Matrix to factorize
    M = k_priv^2 * (inv(V-V*P_val*V)-Vinv)
    # Compute D matrix
    D = dfactor(M; svaltol=1e-6)

    # We return M also, so that the factorization can be redone if needed
    return (D, P_val, X_val, Ω_val, M)
end

"""    (D, P_val, X_val, Ω_val, M) =
            staticInputBlock_DPLQG_ss(Ls, As, Bs, Cs, Vs, Vinvs, Winvs, k_priv, ρ,
                                    m=size(As,1), p=size(Cs,1), r=size(Ls,1))

Compute an input matrix D for the two-block differentially private
steady-state LQG mechanism. D takes linear combinations of
individual signals before privacy-preserving Gaussian noise injection.

Note that we do not guarantee observability of the pair (A,D*C). If
observability is lost and needed, call the dfactor function with the
returned matrix M and a truncation threshold lower than 1e-6 in order
to compute a new matrix D with more rows.

Inputs: one matrix per individual

- Q: size(nusers*m,nusers*m) (m = dim of individual state-space)
- R: size(du,du) (du = dim of control input signal)
- As: size (m,m,nusers)
- Bs: size (m,du,nsusers) - All B matrices pre-multiply the same input u
- Cs: size (p,m,nusers)
- Vs: size (p,p,nusers)
- Vinvs: inverses of Vs (computed externally for modularity/efficiency)
- Winvs: size (m,m,nusers)
- ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.
"""
function staticInputBlock_DPLQG_ss(Q, R, As, Bs, Cs, Ws, Winvs, Vs, Vinvs,
                                ρ, k_priv, nusers=size(As,3), m=size(As,1),
                                p=size(Cs,1), du=size(Bs,2))

    Nm = nusers * m

    A = zeros(Nm, Nm)
    B = zeros(Nm, du)  # the same control input is shared by all agents
    W = zeros(Nm, Nm)

    for i=1:nusers
        A[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = As[:,:,i]
        B[(i-1)*m+1:i*m, :] = Bs[:,:,i]  # stack the Bs in one column
        W[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Ws[:,:,i]
    end

    P = dare(A, B, Q, R)  #  computes ss LQR cost
    N = A'*P*A+Q-P
    L = dfactor(N, svaltol=1e-5)  # factorize L' L = P
    r = size(L,1)

    (D, perf_val, X_val, Ω_val, M) =
        staticInputBlock_DPKF_ss(L, As, Cs, Winvs, Vs, Vinvs, ρ, k_priv, nusers, m, p, r)
    return (D, tr(P*W)+tr(X_val), P, X_val, M)
end



"""    dfactor(M; svaltol=1e-4)

Compute D such that D'D = M for M positive semidefinite, via SVD. Tries to
output a wide D matrix by neglecting the singular values of M smaller than
svaltol*(max svals).
"""
function dfactor(M; svaltol=1e-4)
    (U,S,V) = svd(M)
    svalmax = maximum(S)
    threshold = svaltol*svalmax
    cols = []
    for i=1:length(S)
        if S[i] >= threshold
            push!(cols, i)
        end
    end
    D = zeros(size(M,1), length(cols))
    for j=1:length(cols)
        D[:,j] = U[:,cols[j]] * sqrt(S[cols[j]])
    end
    return sign(D[1,1])*D'
end


"""    evaluateKFperf(D, Ls, As, Cs, Vs, Ws, ρ, k_priv)

Compute the steady-state MSE of the two-block differentially private
Kalman filter, with a static matrix D to combine input signals.
The performance is computed by solving an algebraic Riccati equation.
It corresponds to the MSE at the end of a measurement update state in
the Kalman filter.

Inputs:
- D: size (q,nusers*p)  (q can be arbitrary, p = individual output signal dim.)
one matrix per individual
- Ls: size (r,m,nusers)  (if want estimator of sum_i L_i x_i)
- As: size (m,m,nusers)
- Cs: size (p,m,nusers)
- Ws: size (m,m,nusers) -- process noise cov. matrix
- Vs: size (p,p,nusers) -- measurement noise cov. matrix
- ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.

Outputs: (performance, Σ, Σupdate)
where Σ and Σupdate are the error covariance matrices after the
prediction and measurement update step respectively
"""
function evaluateKFperf(D, Ls, As, Cs, Ws, Vs, ρ, k_priv)
    nusers=size(As,3)
    m=size(As,1)
    p=size(Cs,1)
    r=size(Ls,1)

    Nm = nusers * m
    Np = nusers * p

    L = zeros(r, Nm)
    A = zeros(Nm, Nm)
    C = zeros(Np, Nm)
    W = zeros(Nm, Nm)
    V = zeros(Np, Np)

    for i=1:nusers
        L[:, (i-1)*m+1:i*m] = Ls[:,:,i]
        A[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = As[:,:,i]
        C[(i-1)*p+1:i*p, (i-1)*m+1:i*m] = Cs[:,:,i]
        W[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Ws[:,:,i]
        V[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vs[:,:,i]
        #Winv[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Winvs[:,:,i]
        #Vinv[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vinvs[:,:,i]
    end

    # compute sensitivity
    Δ₂ = 0
    for i=1:nusers
        Δ₂ = max(Δ₂, ρ[i] * maximum(svdvals(D[:,(i-1)*p+1:i*p])))
    end

    V₁ = D*V*D' + k_priv^2  * Δ₂^2 * Matrix{Float64}(I,size(D,1),size(D,1))
    C₁ = D*C
    V₁ = 0.5*(V₁+V₁')  # might have lost perfect symmetry due to numerical issues
    Σ = dare(A', C₁', W, V₁)  #  computes ss cov. after time update step
    Σupdate = Σ - Σ*C₁'*inv(C₁*Σ*C₁'+V₁)*C₁*Σ  # cov. afer meas. update step

    #return (trace(L*Σupdate*L'), Δ₂, Σupdate)
    return tr(L*Σupdate*L'), Σ, Σupdate
end




"""    evaluateLQGperf(D, Q, R, As, Bs, Cs, Ws, Vs, ρ, k_priv)

Compute the steady-state MSE of the two-block differentially private
LQG controller, with a static matrix D to combine input signals.
The performance is computed by solving two algebraic Riccati equations.

Inputs:
- D: size (q,nusers*p)  (q can be arbitrary, p = individual output signal dim.)
one matrix per individual
- Q: size(nusers*m,nusers*m) (m = dim of individual state-space)
- R: size(du,du) (du = dim of control input signal)
- As: size (m,m,nusers)
- Bs: size (m,du,nsusers) - All B matrices pre-multiply the same input u
- Cs: size (p,m,nusers)
- Ws: size (m,m,nusers) -- Process noise
- Vs: size (p,p,nusers) -- Measurement noise
- ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.

Outputs: (performance, Σ, Σupdate, P)
where Σ and Σupdate are the error covariance matrices after the
prediction and measurement update step respectively, P is the solution
of the control algebraic Riccati equation
"""
function evaluateLQGperf(D, Q, R, As, Bs, Cs, Ws, Vs, ρ, k_priv)
    nusers=size(As,3)
    m=size(As,1)
    p=size(Cs,1)
    du=size(R,1)

    Nm = nusers * m
    Np = nusers * p

    A = zeros(Nm, Nm)
    B = zeros(Nm, du)  # the same control input is shared by all agents
    C = zeros(Np, Nm)
    W = zeros(Nm, Nm)
    V = zeros(Np, Np)

    for i=1:nusers
        A[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = As[:,:,i]
        B[(i-1)*m+1:i*m, :] = Bs[:,:,i]  # stack the Bs in one column
        C[(i-1)*p+1:i*p, (i-1)*m+1:i*m] = Cs[:,:,i]
        W[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Ws[:,:,i]
        V[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vs[:,:,i]
    end

    P = dare(A, B, Q, R)  #  computes ss LQR cost matrix

    # compute sensitivity
    Δ₂ = 0
    for i=1:nusers
        Δ₂ = max(Δ₂, ρ[i] * maximum(svdvals(D[:,(i-1)*p+1:i*p])))
    end

    V₁ = D*V*D' + k_priv^2  * Δ₂^2 * Matrix{Float64}(I,size(D,1),size(D,1))
    C₁ = D*C
    V₁ = 0.5*(V₁+V₁')  # might have lost perfect symmetry due to numerical issues
    Σ = dare(A', C₁', W, V₁)  #  computes ss cov. after time update step
    Σupdate = Σ - Σ*C₁'*inv(C₁*Σ*C₁'+V₁)*C₁*Σ  # cov. afer meas. update step

    N = Q+A'*P*A-P
    return tr(P*W + N*Σupdate), P, Σ, Σupdate
end
