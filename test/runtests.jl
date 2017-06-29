#using DynamicDataPrivacy
using PrivateDynamicData
using ControlSystems
using Base.Test

@testset begin
  include("basics_tests.jl")
  include("dpkf_tests.jl")
end
