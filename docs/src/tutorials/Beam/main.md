# Beam Problem

*   An optimal control version of the Singly Supported NonLinear BEAM problem.
*   The energy of a beam of length 1 compressed by a force P is to be minimized.  
*   The control variable is the derivative of the deflection angle.

This problem can be found [here](https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/optcontrol.jl).

## Packages that will be used
```@example Beam
using NLOptControl
nothing # hide
```

## Define and Configure the Problem:
```@example Beam
n=define(numStates=2,numControls=1,XL=[-0.05,-1.0],XU=[-0.05,1.0]);
nothing # hide
```

## Differential Equations
```@example Beam
dx=[:(sin(x2[j])),:(u1[j])]
dynamics!(n,dx)
nothing # hide
```

## Configure the Problem
```@example Beam
configure!(n;(:integrationScheme=>:trapezoidal),(:finalTimeDV=>false),(:tf=>1.0));
nothing # hide
```
## Objective Function
```@example Beam
obj=integrate!(n,:( u1[j]^2 + 350*cos(x2[j]) ) )
@NLobjective(n.ocp.mdl,Min,obj);
nothing # hide
```
## Optimize
```@example Beam
optimize!(n);
nothing # hide
```

## Post Process
```@example Beam
allPlots(n)
```
