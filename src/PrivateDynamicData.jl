module PrivateDynamicData

export gaussianMechConstant, gaussianMechConstant2,
       mean_dp, laplaceMech, gaussianMech,
       staticInputBlock_DPKF_ss, dfactor, evaluateKFperf

include("utils.jl")
include("dpkf.jl")  # Differentially private Kalman filtering

#= Not exported
dpkf_ok: boolean, true if Mosek has been installed (for Kalman filtering)
=#

end # module
