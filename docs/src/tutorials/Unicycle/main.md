# Unicycle Model

### Given:
#### A unicycle trying to get to the goal
#### Dynamic Constraints
$$\dot{x}_1(t)=x_4(t)cos(u_1(t))$$
$$\dot{x}_2(t)=x_4(t)sin(u_1(t))$$
$$\dot{x}_3(t)=u_1(t)$$
$$\dot{x}_4(t)=u_2(t)$$

#### Boundary Conditions
$${x}_1(0)=0 \qquad {x}_1(t_f)=free$$
$${x}_2(0)=pi/2\qquad {x}_2(t_f)=free$$
$${x}_3(0)=0.5\qquad {x}_3(t_f)=free$$
$${x}_4(0)=0\qquad {x}_4(t_f)=free$$

## Find:
#### The control signals that minimize distance to goal $(x_g,y_g)$ within $tplan$
$$J=({x}_1(t_f)-x_g)^2 + ({x}_2(t_f)-y_g)^2)$$

## Solution:
## Packages that will be used
```@example Unicycle
using NLOptControl
nothing # hide
```

## Define the Problem
Next let's write down the boundary conditions into an array:
```@example Unicycle
X0=[0,0,pi/2,0]
XL=[-10,-10,-pi,0]
XU=[10,10,pi,1]
CL=[-1,-3]
CU=[1,3]
nothing # hide
```

## Define the Problem
```@example Unicycle
n = define(numStates=4,numControls=2,X0=X0,XL=XL,XU=XU,CL=CL,CU=CU)
nothing # hide
```

## State and Control Names
```@example Unicycle
names=[:x,:y,:psi,:ux]
descriptions=["X (m)","Y (m)","Yaw Angle (rad)","Longitudinal Velocity (m/s)"]
states!(n,names,descriptions=descriptions)
names = [:r,:ax]
descriptions=["Yaw Rate (rad/s)","Longitudinal Acceleration (m/s^2)"];
controls!(n,names,descriptions=descriptions)
nothing # hide
```

## Differential Equations
```@example Unicycle
dx = [:(ux[j]*cos(psi[j])),:(ux[j]*sin(psi[j])),:(r[j]),:(ax[j])]
dynamics!(n,dx)
nothing # hide
```

## Define and Configure the Problem:
```@example Unicycle
tplan = 7.0
configure!(n;(:Nck=>[50]),(:finalTimeDV=>false), (:tf=>tplan))
nothing # hide
```

## Objective Function
```@example Unicycle
x = n.r.ocp.x[:,1]; y = n.r.ocp.x[:,2]; # pointers to JuMP variables
xg = -2; yg = 4;
@NLobjective(n.ocp.mdl, Min, (x[end]-xg)^2 + (y[end]-yg)^2)
nothing # hide
```

## Optimize
```@example Unicycle
optimize!(n)
nothing # hide
```

## Post Process
```@example Unicycle
plotSettings(;(:size=>(700,700)))
allPlots(n)
```

Taking a closer look at the position:
```@example Unicycle
plotSettings(;(:size=>(400,400)));
statePlot(n,1,1,2;(:lims=>false))
xlims!(-3,2);
ylims!(0,5);
```
