# NLOptControl.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/f480ahs29c85m6ne?svg=true)](https://ci.appveyor.com/project/huckl3b3rry87/nloptcontrol-jl)
[![travis](https://travis-ci.org/JuliaMPC/NLOptControl.jl.svg?branch=master)](https://travis-ci.org/JuliaMPC/NLOptControl.jl)

This software solves **nonlinear control problems** at a **high-level** very **quickly**.

Adds to [juliaOpt](http://www.juliaopt.org/) community by:

* Providing an implementation of direct-collocation methods for solving optimal control problems in julia
* Solving nonlinear optimal control problems at a high-level
* Visualizing the solution

## Documentation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliampc.github.io/NLOptControl.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliampc.github.io/NLOptControl.jl/latest/)

## Installation

If you are using **Linux** make sure that you have **gfortran** to run **Ipopt**:

```bash
sudo apt-get update
sudo apt-get install gfortran
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
```

Also, make sure that you are using at least julia 1.0.0

Then open up julia and install NLOptControl

```julia
Pkg.add("https://github.com/JuliaMPC/NLOptControl.jl")
Pkg.pin("KNITRO",v"0.4")
Pkg.clone("https://github.com/JuliaMPC/NLOptControl.jl")
```

## Citation

If you find [NLOptControl.jl](https://github.com/JuliaMPC/NLOptControl.jl) useful, please cite it:

```LaTeX
@software{nlopt,
  author = {{Huckleberry Febbo}},
  title = {NLOptControl.jl},
  url = {https://github.com/JuliaMPC/NLOptControl.jl},
  version = {0.0.1},
  date = {2017-06-17},
}
```

## Acknowledgements

* [JuMP.jl](https://jump.readthedocs.io/en/latest/) is an important part of this NLOptControl.jl and discussions with Miles Lubin where helpful
* Chris Rackauckas is a very helpful member of the julia community and has provided me support and advice multiple times his software [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) is also part of NLOptControl.jl
