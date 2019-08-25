PrivateDynamicData.jl
=====================

*Differentially private dynamic data analysis in Julia*

---

## Set-up ##

This package is intended to work with Julia v1.0 and later.
To use with Julia prior to v1.0, please see the branch julia0.5 of this repository, as a starting point.

This is branch jump0.18 of the repository: it is intended to work with [JuMP](https://github.com/JuliaOpt/JuMP.jl) version 0.18 (branch release-0.18), i.e., prior to the change from MathProgBase to MathOptInterface. It also uses [Mosek](https://www.mosek.com/) version 8.1 to solve semidefinite programs, and the corresponding branch b0.9 of the [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl) package. These correct version of these packages are checked-out automatically by the package manager, but you need to have Mosek installed and a valid license. See the instruction for the [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl) package. You can also refer to the JuMP documentation on
[installing solvers](http://www.juliaopt.org/JuMP.jl/v0.18/installation.html#getting-solvers).

Note that Mosek requires a license, but free Academic licenses are available.
If Mosek is not installed, this package can still be installed but the functions for differentially private Kalman filtering will not work. Adding support for other semidefinite solvers than Mosek or other versions of Mosek and streamlining the installation process is currently WIP in other branches of this package.

To add the package, in the package manager (press ] in the REPL)
```julia
pkg> add "https://github.com/cipherlab-poly/PrivateDynamicData.jl.git"#jump0.18
```

To use the package functionalities
```julia
using PrivateDynamicData
```

## Run Tests ##

```julia
pkg> test PrivateDynamicData
```

## Functions ##

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

## Contributors ##

* Jerome Le Ny
