@testset "Differential Private Kalman Filtering Tests" begin
  @testset "Static Input Block Design" begin
    ϵ = log(3); δ = 0.05
    k_priv = gaussianMechConstant(ϵ,δ)

    nusers = 10; ρ = ones(1,nusers)
    Ls = ones(1,1,nusers); As = 0.9*ones(1,1,nusers); Cs = ones(1,1,nusers)
    Vs = 0.01*ones(1,1,nusers); Ws = 0.2*ones(1,1,nusers)
    Vinvs = zero(Vs); Winvs = zero(Ws)
    for i=1:nusers
    	Vinvs[:,:,i] = inv(Vs[:,:,i])
    	Winvs[:,:,i] = inv(Ws[:,:,i])
    end
    (D, P_val, X_val, Ω_val, M) =
      staticInputBlock_DPKF_ss(Ls, As, Cs, Winvs, Vs, Vinvs, ρ, k_priv)
    cost = evaluateKFperf(D, Ls, As, Cs, Ws, Vs, ρ, k_priv)[1]

    @test sum(D) ≈ nusers atol=1e-6 # D should be all ones
    @test cost ≈ 1.6244783332772452
    @test tr(X_val) ≈ cost atol=1e-6

    # LQG test - need to finish, writing evaluation function for LQG
    Bs = ones(1,1,nusers)
    Q = 1.0 * Matrix(I,nusers,nusers); R = 0.5
    (Dlqg, perflqg, Plqg, X_val_lqg, Mlqg) =
      staticInputBlock_DPLQG_ss(Q, R, As, Bs, Cs, Ws, Winvs, Vs, Vinvs, ρ, k_priv)
    costlqg = evaluateLQGperf(Dlqg, Q, R, As, Bs, Cs, Ws, Vs, ρ, k_priv)[1]
    @test perflqg ≈ costlqg atol=1e-6
  end
end
