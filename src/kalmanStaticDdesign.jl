using JuMP
using Mosek

"""    computeDsteadyState(L::Matrix, A::Matrix, C::Matrix, V::Matrix,
                           Winv::Matrix, α::Real, ρ::Real, Vinv=inv(V),
                           m=size(A,1), p=size(C,1), r=size(L,1))
Compute a static D matrix for signal preparation of one user, to be
used for differentially private Kalman filtering.
"""
function computeDsteadyState(L, A, C, V, Winv, α, ρ, Vinv=inv(V),
                             m=size(A,1), p=size(C,1), r=size(L,1))
  modl = Model(solver=MosekSolver())

  @variable(modl, X[1:r, 1:r], SDP)
  @variable(modl, Ω[1:m, 1:m], SDP)
  @variable(modl, P[1:p, 1:p], SDP)

  @objective(modl, Min, trace(X))

  @SDconstraint(modl, hvcat((2,2), X, L, L', Ω) >= zeros(r+m,r+m))
  pibound = inv(α^2 * eye(p) + V);
  @SDconstraint(modl, P <= pibound)

  @SDconstraint(modl, hvcat((2,2),
    C'*P*C-Ω+Winv, Winv*A, A'*Winv, Ω+A'*Winv*A) >= zeros(2m,2m))

  status = solve(modl)

  #println("Solution status: ", status)
  #println("Objective value: ", getobjectivevalue(modl))
  P_val = getvalue(P)
  #=
  X_val = getvalue(X)
  Ω_val = getvalue(Ω)
  Σ_val = inv(Ω_val)
  println("Difference between Σ^{-1} and (A Σ A^T + W)^{-1}+P")
  E2 = inv(Σ_val)-inv(A*Σ_val*A'+W)-P_val
  norm(E2)
  =#

  # Matrix to factorize
  M₁ = α^2 * (inv(V-V*P_val*V)-Vinv)
  return chol(M₁)/ρ
  #return D'/ρ
end
