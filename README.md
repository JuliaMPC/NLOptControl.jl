# NLOptControl.jl


[![Travis](https://travis-ci.org/JuliaMPC/NLOptControl.jl.svg?branch=master)](https://travis-ci.org/JuliaMPC/NLOptControl.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliampc.github.io/MPCDocs.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliampc.github.io/MPCDocs.jl/latest/)


A current limitation of optimization modeling software, such as [JuMP](https://github.com/JuliaOpt/JuMP.jl) is that it does not allow for ease of adding integral constraint equations. The aim of this package is to provide functionality to:
  * Allow the user to write a system of coupled ODEs that governs the dynamic behavior of the model

  * provide an implementation of the [hp-pseudospectral method](http://vdol.mae.ufl.edu/JournalPublications/TOMS-GPOPS-II-August-2013.pdf) written in written in the [Julia-language](http://julialang.org/)

## Documentation

The full documentation can be found [here](https://juliampc.github.io/MPCDocs.jl/stable/).

## Installation

In Julia, you can install the NLOptControl.jl package by typing:
```julia
Pkg.clone("https://github.com/JuliaMPC/PrettyPlots.jl")
Pkg.clone("https://github.com/JuliaMPC/VehicleModels.jl")
Pkg.clone("https://github.com/JuliaMPC/NLOptControl.jl")
```

## Tutorial

-> currently out of date!
All Examples can be:

  *  Viewed remotely on  using the [jupyter nbviewer](http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/).
  *  Viewed locally and interacted using IJulia

      To do this in julia type:
      ```julia
      using IJulia
      notebook(dir=Pkg.dir("NLOptControl/examples/"))
      ```

## Citation

If you find this package useful, please cite this paper (accepted but not yet released to public)

```
@Conference{Febbo2017,
  author    = {Huckleberry Febbo, Jiechao Liu, Paramsothy Jayakumar, Jeffrey L. Stein, Tulga Ersal},
  title     = {Moving Obstacle Avoidance for Large, High-Speed Autonomous Ground Vehicles},
  year      = {2017},
  publisher = {IEEE},
}
```
