# Bryson Denham

This problem can be found [here](http://www.gpops2.com/Examples/Bryson-Denham.html).

## Packages that will be used
```@example BrysonDenham
using NLOptControl
nothing # hide
```

## Define the Problem:
```@example BrysonDenham
n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[1/9,NaN]);
nothing # hide
```

## Differential Equations
```@example BrysonDenham
dx=[:(x2[j]),:(u1[j])]
dynamics!(n,dx)
nothing # hide
```

## Configure the Problem
```@example BrysonDenham
configure!(n;(:Nck=>[100]),(:finalTimeDV=>true));
nothing # hide
```

## Objective Function
```@example BrysonDenham
obj=integrate!(n,:(0.5*u1[j]^2));
@NLobjective(n.ocp.mdl,Min,obj);
nothing # hide
```
## Optimize
```@example BrysonDenham
optimize!(n);
nothing # hide
```

## Post Process
```@example BrysonDenham
allPlots(n)
```
