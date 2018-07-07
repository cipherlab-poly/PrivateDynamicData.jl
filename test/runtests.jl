#using DynamicDataPrivacy
using PrivateDynamicData
using ControlSystems
using Base.Test

@testset begin
  include("basics_tests.jl")
  if dpkf_ok
      include("dpkf_tests.jl")
  end
end
