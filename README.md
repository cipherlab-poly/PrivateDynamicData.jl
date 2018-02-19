PrivateDynamicData.jl
=====================

*Differentially private dynamic data analysis in Julia*

---

## Set-up ##

```julia
Pkg.clone("https://github.com/cipherlab-poly/PrivateDynamicData.jl.git")
using PrivateDynamicData
```

Currently the package requires Mosek and Mosek.jl, to solve semidefinite programs.
Simplifying the installation process is WIP.

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
