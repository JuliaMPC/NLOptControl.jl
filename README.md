# NLOptControl.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/f480ahs29c85m6ne?svg=true)](https://ci.appveyor.com/project/huckl3b3rry87/nloptcontrol-jl)
[![travis](https://travis-ci.org/JuliaMPC/NLOptControl.jl.svg?branch=master)](https://travis-ci.org/JuliaMPC/NLOptControl.jl)
[![NLOptControl](http://pkg.julialang.org/badges/NLOptControl_0.6.svg)](http://pkg.julialang.org/detail/NLOptControl)
[![NLOptControl](http://pkg.julialang.org/badges/NLOptControl_0.7.svg)](http://pkg.julialang.org/detail/NLOptControl)

This software solves **nonlinear control problems** at a **high-level** very **quickly**.

Adds to [juliaOpt](http://www.juliaopt.org/) community by:
 * Providing an implementation of the direct-collocation methods for solving optimal control problems in julia
 * Solving nonlinear optimal control problems at a high-level
 * Visualizing the solution

## Documentation
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliampc.github.io/NLOptControl.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliampc.github.io/NLOptControl.jl/latest/)

## Installation
```julia
Pkg.add("NLOptControl")
```

Either [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) or [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl) can be used for the nonlinear solver. Since **Ipopt** is free, all of the examples will use it and it can be added with:
```julia
Pkg.add("Ipopt");Pkg.build("Ipopt");
```

If you are using **Linux** make sure that you have **gfortran** to run **Ipopt**:
```
$sudo apt-get update
$sudo apt-get install gfortran
$sudo apt-get liblapack-dev
$sudo apt-get libblas-dev
```

## Add Ons
For additional functionality, there a couple other packages that can be installed:
```julia
Pkg.add("VehicleModels")
Pkg.add("PrettyPlots")
Pkg.clone("https://github.com/JuliaMPC/MichiganAutonomousVehicles.jl")
```

Also, a plotting backend will be required and there are several options:
```julia
Pkg.add("PGFPlots");Pkg.build("PGFPlots")   # best looking
Pkg.add("GR");Pkg.build("GR");              # most reliable
Pkg.add("PyPlot");Pkg.build("PyPlot")       # also a good option  
```

## Citation

If you find [NLOptControl.jl](https://github.com/JuliaMPC/NLOptControl.jl) useful, please cite it:
```
@software{nlopt,
  author = {{Huckleberry Febbo}},
  title = {NLOptControl.jl},
  url = {https://github.com/JuliaMPC/NLOptControl.jl},
  version = {0.0.1},
  date = {2017-06-17},
}
```

If you find [VehicleModels.jl](https://github.com/JuliaMPC/VehicleModels.jl) useful, please cite this paper:
```
@Conference{Febbo2017,
  author    = {Huckleberry Febbo, Jiechao Liu, Paramsothy Jayakumar, Jeffrey L. Stein, Tulga Ersal},
  title     = {Moving Obstacle Avoidance for Large, High-Speed Autonomous Ground Vehicles},
  year      = {2017},
  publisher = {IEEE}
}
```

## Acknowledgements
* [JuMP.jl](https://jump.readthedocs.io/en/latest/) is an important part of this NLOptControl.jl and discussions with Miles Lubin where helpful
* Chris Rackauckas is a very helpful member of the julia community and has provided me support and advice multiple times his software [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) is also part of NLOptControl.jl
