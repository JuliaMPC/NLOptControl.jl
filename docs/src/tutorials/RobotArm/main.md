# RobotArm

This problem can be found [here](http://www.gpops2.com/Examples/RobotArm.html).
 although that example is missing initial and final state constraints and limits on x4


## Packages that will be used
```@example RobotArm
using NLOptControl
nothing # hide
```

## Define the Problem:
```@example RobotArm
n = define(numStates=6,numControls=3,X0=[9/2,0.0,0.0,0.0,pi/4,0.0],XF=[9/2,0.0,2*pi/3,0.0,pi/4,0.0],XL=[NaN,NaN,NaN,0.0,NaN,NaN],XU=[NaN,NaN,NaN,1.0,NaN,NaN],CL=[-1.,-1.,-1.],CU=[1.,1.,1.])
nothing # hide
```

## Constants
```@example RobotArm
EP = 2*eps(); # to avoid divide/0
Q = 5
nothing # hide
```

# Differential Equations

```@example RobotArm
# expressions
I_t = :((($Q-x1[j])^3+x1[j]^3)/3*sin(x5[j])^2)
I_p = :((($Q-x1[j])^3+x1[j]^3)/3 )

# Diff Eqs
dx = Array{Expr}(6,)
dx[1] = :(x2[j])
dx[2] = :(u1[j]/$Q)
dx[3] = :(x4[j])
dx[4] = :(u2[j]/($I_t+$EP))
dx[5] = :(x6[j])
dx[6] = :(u3[j]/($I_p+$EP))
```
Then add the differential equations to the model:
```@example RobotArm
dynamics!(n,dx)
```

## Configure the Problem:
```@example RobotArm
configure!(n;(:finalTimeDV=>true))
nothing # hide
```

## Objective Function
```@example RobotArm
@NLobjective(n.ocp.mdl,Min,n.ocp.tf)
nothing # hide
```

## Optimize
```@example RobotArm
optimize!(n)
nothing # hide
```

## Post Process
```@example RobotArm
allPlots(n)
```
