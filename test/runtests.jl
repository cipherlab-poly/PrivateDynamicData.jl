#  Copyright 2019, Jerome Le Ny

using PrivateDynamicData

using Test  # Note: test-specific dependency

using Random
using Distributions

using ControlSystems
using LinearAlgebra

@testset "Complete Test Set" begin
  include("basics_tests.jl")
  if PrivateDynamicData.dpkf_ok
      include("dpkf_tests.jl")
  end
end
