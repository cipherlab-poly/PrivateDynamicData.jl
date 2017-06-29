#__precompile__()

module PrivateDynamicData

export gaussianMechConstant, gaussianMechConstant2,
       laplaceMech, gaussianMech,
       staticInputBlock_DPKF_ss, dfactor, evaluateKFperf

include("utils.jl")
include("dpkf.jl")  # Differentially private Kalman filtering

end  # module
