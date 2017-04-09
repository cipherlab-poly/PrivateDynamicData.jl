# test.jl
# Solves the matrix design problem for one agent in order to
# enforce differential privacy in Kalman filtering

using JuMP
using Mosek

# Problem data
Ts = 1;
A = [1 Ts;0 1];
# A = [0 1;0 1];

n = size(A,1);  # state space dimension

# B = [Ts^2/2;Ts];
# #sig = 1;
#% #W = B*sig*B';
# #W = eye(2);
# W = 0.5*randn(2); W = W*W';
W = 1e-1 * eye(n);

C = eye(n);
m = size(C,1);

Winv = inv(W);
V = eye(m);

# Privacy multiplicative constant
alpha = 1.7;

# Solve the SDP
pibound = inv(alpha^2 * eye(m) + V);

mod = Model()
@variable(mod, X[1:n, 1:n], SDP)
@variable(mod, Ω[1:n, 1:n], SDP)
@variable(mod, P[1:m, 1:m], SDP)

@objective(mod, Min, trace(X))

@SDconstraint(mod, hvcat((2,2), X, eye(n), eye(n), Ω) >= zeros(2n,2n))
@SDconstraint(mod, P <= pibound)
@SDconstraint(mod, hvcat((2,2),
  C'*P*C-Ω+Winv, Winv*A, A'*Winv, Ω+A'*Winv*A) >= zeros(2n,2n))

# Solver, recover and analyze the solution
status = solve(mod)

println("Solution status: ", status)
println("Objective value: ", getobjectivevalue(mod))

X_val = getvalue(X)
Ω_val = getvalue(Ω)
P_val = getvalue(P)

Σ_val = inv(Ω_val)

println("Difference between Σ^{-1} and (A Σ A^T + W)^{-1}+P")
E2 = inv(Σ_val)-inv(A*Σ_val*A'+W)-P_val
norm(E2)

println("Matrix to factorize")
M1 = alpha^2 * (inv(V-V*P_val*V)-inv(V))
