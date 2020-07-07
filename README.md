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
```
sudo apt-get update
sudo apt-get install gfortran
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
```

Also, make sure that you are using julia 0.6.4
```
sudo wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz
sudo tar -xvf julia-0.6.4-linux-x86_64.tar.gz -C /opt
```

Then open up julia and install a few specific packages and NLOptControl
```julia
Pkg.add("ReverseDiffSparse")
Pkg.checkout("ReverseDiffSparse")
Pkg.add("KNITRO")
Pkg.pin("KNITRO",v"0.4")
Pkg.clone("https://github.com/JuliaMPC/NLOptControl.jl")
```

## Citation
If you find [NLOptControl.jl](https://github.com/JuliaMPC/NLOptControl.jl) useful, please cite it:
```
@misc{febbo2020nloptcontrol,
    title={NLOptControl: A modeling language for solving optimal control problems},
    author={Huckleberry Febbo and Paramsothy Jayakumar and Jeffrey L. Stein and Tulga Ersal},
    year={2020},
    eprint={2003.00142},
    archivePrefix={arXiv},
    primaryClass={cs.MS}
}
```

## Acknowledgements
* [JuMP.jl](https://jump.readthedocs.io/en/latest/) is an important part of this NLOptControl.jl and discussions with Miles Lubin where helpful
* Chris Rackauckas is a very helpful member of the julia community and has provided me support and advice multiple times his software [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) is also part of NLOptControl.jl
