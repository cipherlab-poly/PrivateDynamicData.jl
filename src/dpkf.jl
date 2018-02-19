using JuMP
using Mosek
#using SCS
using ControlSystems


"""
    (D, P_val, X_val, Ω_val, M) =
        staticInputBlock_DPKF_ss(Ls, As, Cs, Vs, Vinvs, Winvs, k_priv, ρ,
                                 m=size(As,1), p=size(Cs,1), r=size(Ls,1))

Compute an input matrix D for the two-block differentially private
steady-state Kalman filter mechanism. D takes linear combinations of
individual signals before privacy-preserving Gaussian noise injection.

Inputs: one matrix per individual

- Ls: size (r,m,nusers)  (if want estimator of sum_i L_i x_i)
- As: size (m,m,nusers)
- Cs: size (p,m,nusers)
- Vs: size (p,p,nusers)
-	Vinvs: inverses of Vs (computed externally for modularity/efficiency)
-	Winvs: size (m,m,nusers)
-	ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.
"""
function staticInputBlock_DPKF_ss(Ls, As, Cs, Vs, Vinvs, Winvs, ρ, k_priv,
								  nusers=size(As,3), m=size(As,1),
                                  p=size(Cs,1), r=size(Ls,1))
    modl = Model(solver=MosekSolver())
    #modl = Model(solver=SCSSolver())

    Nm = nusers * m
    Np = nusers * p

    L = zeros(r, Nm)
    A = zeros(Nm, Nm)
    C = zeros(Np, Nm)
    V = zeros(Np, Np)
    Winv = zeros(Nm, Nm)
    Vinv = zeros(Np, Np)

    for i=1:nusers
      L[:, (i-1)*m+1:i*m] = Ls[:,:,i]
	  A[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = As[:,:,i]
	  C[(i-1)*p+1:i*p, (i-1)*m+1:i*m] = Cs[:,:,i]
	  V[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vs[:,:,i]
	  Winv[(i-1)*m+1:i*m, (i-1)*m+1:i*m] = Winvs[:,:,i]
	  Vinv[(i-1)*p+1:i*p, (i-1)*p+1:i*p] = Vinvs[:,:,i]
    end

    @variable(modl, X[1:r, 1:r], SDP)
    @variable(modl, Ω[1:Nm, 1:Nm], SDP)
    @variable(modl, P[1:Np, 1:Np], SDP)

    @objective(modl, Min, trace(X))

    @SDconstraint(modl, hvcat((2,2), X, L, L', Ω) >= zeros(r+Nm,r+Nm))

    @SDconstraint(modl, hvcat((2,2),
    C'*P*C-Ω+Winv, Winv*A, A'*Winv, Ω+A'*Winv*A) >= zeros(2*Nm,2*Nm))

    for i=1:nusers
	    e = zeros(p,Np); e[:,(i-1)*p+1:i*p]=eye(p)
	    @SDconstraint(modl, hvcat((2,2),
	    eye(p)/(k_priv^2*ρ[i]^2)+Vinvs[:,:,i], e, e', V-V*P*V) >= zeros(Np+p,Np+p))
    end

    status = solve(modl)

    println("============================================")
    println("Solution status: ", status)
    println("Objective value: ", getobjectivevalue(modl))
    println("============================================")

    P_val = getvalue(P)

    X_val = getvalue(X)
    Ω_val = getvalue(Ω)

    ## Debug
    #Σ_val = inv(Ω_val)
    #println("Norm of the difference between Σ^{-1} and (A Σ A' + W)^{-1}+C' P C")
    #E1 = C'*P_val*C+inv(A*Σ_val*A'+inv(Winv))-inv(Σ_val)
    #println(norm(E1))
    ##

    # Matrix to factorize
    M = k_priv^2 * (inv(V-V*P_val*V)-Vinv)
    # Compute D matrix
    D = dfactor(M; svaltol=1e-3)

    #return chol(Hermitian(M₁))
    #return (M, P_val, Ω_val, Σ_val, X_val)
    #	return (M, P_val, Ω_val, Σ_val, X_val, A, C, Winv, V)
    return (D, P_val, X_val, Ω_val, M)
end


"""
    dfactor(M; svaltol=1e-3)

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


"""
    evaluateKFperf(D, Ls, As, Cs, Vs, Ws, ρ, k_priv)

Compute the steady-state MSE of the two-block differentially private
Kalman filter, with a static matrix D to combine input signals.
The performance is computed by solving an algebraic Riccati equation.
It corresponds to the MSE at the end of a measurement update state in
the Kalman filter.

Inputs: one matrix per individual

- D: size (q,nusers*p)  (q can be arbitrary, p = individual output signal dim.)
- Ls: size (r,m,nusers)  (if want estimator of sum_i L_i x_i)
- As: size (m,m,nusers)
- Cs: size (p,m,nusers)
- Vs: size (p,p,nusers)
-	Ws: size (m,m,nusers)
-	ρ: vector of size nusers, for the adjacency definition: ||yᵢ-yᵢ'||₂ <= ρᵢ
for user i.
- k_priv: proportionality constant for Gaussian privacy-preserving noise
injection. Compute it with k_priv = gaussianMechConstant(ϵ,δ) for your
choice of ϵ, δ.
"""
function evaluateKFperf(D, Ls, As, Cs, Vs, Ws, ρ, k_priv)
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

	V₁ = D*V*D' + k_priv^2  * Δ₂^2 * eye(size(D,1))
	C₁ = D*C
	Σ = dare(A', C₁', W, V₁)  #  computes ss cov. after time update step
	Σupdate = Σ - Σ*C₁'*inv(C₁*Σ*C₁'+V₁)*C₁*Σ  # cov. afer meas. update step

    #return (trace(L*Σupdate*L'), Δ₂, Σupdate)
	return trace(L*Σupdate*L')
end
