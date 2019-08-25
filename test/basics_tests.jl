@testset "Basic Utilities Tests" begin
  @testset "Utils" begin
    ϵ = log(3); δ = 0.05

    @test gaussianMechConstant(ϵ, δ) ≈ 1.7563398731147595 atol=1e-16
    tmp = gaussianMechConstant2(ϵ, δ)
    @test tmp[1] ≈ 1.7563398731147595 atol=1e-16
    @test tmp[2] ≈ 2.3095249420840473 atol=1e-16
  end

  @testset "Basic mechanisms" begin
    ϵ = log(2); δ = 0.04

    Random.seed!(12345)  # reseed the RNG
    D = Random.rand(100)
    @test gaussianMech(D, Distributions.mean, 2.0/100, ϵ, δ)[1] == 0.4764940748202759
    @test laplaceMech(D, Distributions.mean, 2.0/100, ϵ)[1] == 0.5366973705020675
    @test mean_dp(D, ϵ, δ) == 0.47514367413387
    @test mean_dp(D, ϵ, 0.0) == 0.4920613726730682
    @test mean_dp(D, ϵ) == 0.4886977015049402
  end
end
