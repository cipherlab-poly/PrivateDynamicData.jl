module PrivateDynamicData

export gaussianMechConstant, gaussianMechConstant2,
       laplaceMech, gaussianMech,
       dpkf_ok,  # constant to check if the SDP solver for Kalman filtering has been installed
       staticInputBlock_DPKF_ss, dfactor, evaluateKFperf

include("utils.jl")
include("dpkf.jl")  # Differentially private Kalman filtering

end # module
