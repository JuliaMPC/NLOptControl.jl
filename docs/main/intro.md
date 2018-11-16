# NLOptControl.jl

This software solves **nonlinear control problems** at a **high-level** very **quickly**.

Adds to [juliaOpt](http://www.juliaopt.org/) community by:
 * Providing an implementation of direct-collocation methods for solving optimal control problems in julia
 * Solving nonlinear optimal control problems at a high-level
 * Visualizing the solution

## Installation
```julia
Pkg.add("NLOptControl")
```

If you are using **Linux** make sure that you have **gfortran** to run **Ipopt**:
```
$sudo apt-get update
$sudo apt-get install gfortran
$sudo apt-get liblapack-dev
$sudo apt-get libblas-dev
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
