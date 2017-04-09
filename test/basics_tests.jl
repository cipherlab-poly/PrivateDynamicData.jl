@testset "Differential Privacy Tests" begin
  @testset "Basics" begin
    ϵ = log(3); δ = 0.05

    @test gaussianMechConstant(ϵ, δ) == 1.756339873114760
    @test gaussianMechConstant2(ϵ, δ) == (1.756339873114760, 2.309524942084048)
  end

  @testset "Basic mechanisms" begin
    ϵ = log(2); δ = 0.01

    @test 1 == 1
  end
end
