PrivateDynamicData.jl
=====================

*Differentially private dynamic data analysis in Julia*

---

## Set-up ##

This package is intended to work with Julia v1.0 and later.
To use with Julia prior to v1.0, please see the version tagged 0.0.1 (initially developed for Julia 0.5). The functions in this package have been (very lightly) tested with Julia v1.3.

To add the package, in the package manager (press ] in the REPL)
```julia
pkg> add "https://github.com/cipherlab-poly/PrivateDynamicData.jl.git"
```

To use the package functionalities
```julia
using PrivateDynamicData
```

The functions for differentially private Kalman filtering require solving semidefinite programs. For simplicity, we adopt a rigid procedure for the installation of specific versions of the necessary optimization packages.  This package has been tested with the [Mosek](https://www.mosek.com/) solver, version 9.1.10. Installing Mosek is necessary to use this version of the package and it is recommended to try installing this specific version first. Note that Mosek requires a license, but free Academic licenses are available. If Mosek is not installed, this package can still be installed but the corresponding functions will not work. Adding support for other semidefinite programming solvers than Mosek should be fairly straightforward but is left to the user for now.

This package will automatically install [JuMP](http://www.juliaopt.org/JuMP.jl/v0.19.2/), version 0.20.1 (this version uses the [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)). Since we rely on Mosek to solve the semidefinite programs, it will also install [MosekTools.jl](https://github.com/JuliaOpt/MosekTools.jl), version 0.9.1. If the installation does not work properly, you might try to work in a clean Julia [environment](https://julialang.github.io/Pkg.jl/v1/environments/), or refer to the JuMP documentation on
[installing solvers](http://www.juliaopt.org/JuMP.jl/v0.19.2/installation/) and
also [MosekTools](https://github.com/JuliaOpt/MosekTools.jl).

## Run Tests ##

```julia
pkg> test PrivateDynamicData
```

## Functions Available ##

Use the help to access the documentation of the following functions.

Utilities:
* gaussianMechConstant
* gaussianMechConstant2
* mean_dp
* laplaceMech
* gaussianMech

Differentially Private Kalman Filtering:
* staticInputBlock_DPKF_ss
* dfactor
* evaluateKFperf

## Examples ##

A few examples showing how to use the functions in this package are available in the form of a Jupyter notebook [here](https://github.com/jleny/DifferentialPrivacy-course).

## Contributors ##

* Jerome Le Ny
