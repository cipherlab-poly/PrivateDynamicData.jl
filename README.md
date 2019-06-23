PrivateDynamicData.jl
=====================

*Differentially private dynamic data analysis in Julia*

---

## Set-up ##

This package is intended to work with Julia v1.0 and later.
To use with Julia prior to v1.0, please see the other branches of this repository.

To add the package, in the package manager (press ] in the REPL)
```julia
pkg> add "https://github.com/cipherlab-poly/PrivateDynamicData.jl.git"
```

To use the package functionalities
```julia
using PrivateDynamicData
```

Currently the functions in the package for differentially private Kalman filtering 
require [Mosek](https://www.mosek.com/) interfaced with 
[JuMP](http://www.juliaopt.org/JuMP.jl/v0.19.2/), vs. 0.19 or later (i.e., using
the [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)).
Mosek is used here to solve semidefinite programs. 
As the installation instructions for the Julia-Mosek interface tend to change over time 
and depend on the version of Mosek that you are using, you might have to 
adapt the current installation process for the
[Mosek.jl](https://github.com/JuliaOpt/Mosek.jl) and
[MosekTools.jl](https://github.com/JuliaOpt/MosekTools.jl)
packages to correspond to your version of Mosek.
This package currently uses Mosek 0.9 and hence the master branch of these interfaces
```julia
pkg> add Mosek#master
pkg> add MosekTools#master
```
Please refer to the JuMP documentation on 
[installing solvers](http://www.juliaopt.org/JuMP.jl/v0.19.2/installation/) and 
also [MosekTools](https://github.com/JuliaOpt/MosekTools.jl).

Note that Mosek requires a license, but free Academic licenses are available. 
If Mosek is not installed, this package can still be installed but the corresponding 
functions will not work. Adding support for other semidefinite solvers than Mosek and
streamlining the installation process is currently WIP.

## Run Tests ##

```julia
pkg> test PrivateDynamicData
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
