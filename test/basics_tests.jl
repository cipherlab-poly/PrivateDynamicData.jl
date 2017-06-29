@testset "Basic Utilities Tests" begin
  @testset "Utils" begin
    ϵ = log(3); δ = 0.05

    @test gaussianMechConstant(ϵ, δ) == 1.756339873114760
    @test gaussianMechConstant2(ϵ, δ) == (1.756339873114760, 2.309524942084048)
  end

  @testset "Basic mechanisms" begin
    ϵ = log(2); δ = 0.04

    srand(12345)  # reseed the RNG
    D = rand(100)
    @test gaussianMech(D, mean, 2.0/100, ϵ, δ)[1] == 0.4764940748202759
    @test laplaceMech(D, mean, 2.0/100, ϵ)[1] == 0.5366973705020675
  end
end
