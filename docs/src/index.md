# NLOptControl.jl Documentation


## Introduction

This software solves **nonlinear control problems** at a **high-level** very **quickly**.

It added to [juliaOpt](http://www.juliaopt.org/) community by:
 * Providing an implementation of of the [hp-pseudospectral method](http://vdol.mae.ufl.edu/JournalPublications/TOMS-GPOPS-II-August-2013.pdf) written in julia
 * Incorporating model predictive control functionality
 * Automatically visualizing the solution

## Installation

To install
```julia
Pkg.add("NLOptControl")
```

If you are using **Linux** make sure that you have **gfortran** to run **Ipopt**:
```julia
sudo apt-get install gfortran liblapack-dev libblas-dev
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

## 2017 juliaCon Workshop Notebook (OUT OF DATE!)

After installation, the notebook can be viewed:
```julia
using IJulia
notebook(dir=Pkg.dir("NLOptControl/examples"))
```

Also, on the left side of this site, there are many tutorials that provide complete examples for using this software. Please look at these for information on how to use this tool.


## Acknowledgements
* [JuMP.jl](https://jump.readthedocs.io/en/latest/) is an important part of this NLOptControl.jl and discussions with Miles Lubin where helpful
* Chris Rackauckas is a very helpful member of the julia community and has provided me support and advice multiple times his software [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) is also part of NLOptControl.jl

## Exported Functions

The following link provides documentation all of the exported functions for `NLOptControl.jl`

```@contents
Pages=[
    "functions/NLOptControl.md"
    ]
Depth=1
```
