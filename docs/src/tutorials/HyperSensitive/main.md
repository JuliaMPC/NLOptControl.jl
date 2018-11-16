# HyperSensitive

This problem can be found [here](http://www.gpops2.com/Examples/Brachistochrone.html).


## Packages that will be used
```@example HyperSensitive
using NLOptControl
nothing # hide
```

## Define the Problem:
```@example HyperSensitive
n=define(numStates=1,numControls=1,X0=[1.5],XF=[1.])
nothing # hide
```
## Differential Equations
```@example HyperSensitive
dx=[:(-x1[j]^3+u1[j])]
dynamics!(n,dx)
nothing # hide
```

## Configure the Problem:
```@example HyperSensitive
configure!(n;(:Nck=>[3,3,3,3,3,3,3,3,3,3,3,3]),(:finalTimeDV=>false),(:tf=>10000.0))
nothing # hide
```

## Objective Function
```@example HyperSensitive
obj=integrate!(n,:( 0.5*x1[j]^2 + 0.5*u1[j]^2) )
@NLobjective(n.ocp.mdl,Min,obj);
```

## Optimize
```@example HyperSensitive
optimize!(n);
nothing # hide
```
## Post Process
```@example HyperSensitive
plotSettings(;(:size=>(1200,1200)));
allPlots(n)
```
