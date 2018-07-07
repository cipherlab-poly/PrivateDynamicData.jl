PrivateDynamicData.jl
=====================

*Differentially private dynamic data analysis in Julia*

---

## Set-up ##

```julia
Pkg.clone("https://github.com/cipherlab-poly/PrivateDynamicData.jl.git")
using PrivateDynamicData
```

Currently the functions in the package for differentially private Kalman filtering require [Mosek](https://www.mosek.com/) and [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl), to solve semidefinite programs. Mosek.jl should be installed PRIOR to installing this package, using

```julia
Pkg.add("Mosek")
```

Note that Mosek requires a license (free Academic licenses are available). If Mosek.jl is not installed, this package can still be installed but the corresponding functions will not work. Replacing Mosek by a free, open-source alternative is currently WIP.

## Run Tests ##

```julia
Pkg.test("PrivateDynamicData")
```

## Functions ##

Use the help to access the documentation of the following functions.

Utilities:
* gaussianMechConstant
* gaussianMechConstant2
* laplaceMech
* gaussianMech

Differentially Private Kalman Filtering:
* staticInputBlock_DPKF_ss
* dfactor
* evaluateKFperf

## Contributors ##

* Jerome Le Ny
