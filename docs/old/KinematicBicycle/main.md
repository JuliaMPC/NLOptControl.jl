# Kinematic Bicycle Model

The vehicle model comes from the [BARC-project](https://github.com/MPC-Berkeley/barc)

## Packages that will be used
```@example Bicycle
using NLOptControl, Parameters
nothing # hide
```

## Parameters form VehicleModels.jl
```@example Bicycle
pa = Vpara(x0_=0.)  
@unpack_Vpara pa # vehicle parameters
X0 = [x0_,y0_,psi0_,u0_]
XF = [NaN,NaN,NaN,NaN]
XL = [x_min,y_min,psi_min,u_min]
XU = [x_max,y_max,psi_max,u_max]
CL = [sa_min,ax_min]
CU = [sa_max,ax_max]
nothing # hide
```

## Define the Problem
```@example Bicycle
n = define(numStates=4,numControls=2,X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU)
n.ocp.params = [pa] # add vehicle parameters
nothing # hide
```

## Define Case
```@example Bicycle
c = Dict("goal"=>Dict("x"=>0.,"yVal"=>100.,"tol"=>5.), "obstacle"=>Dict("x0"=>[0.],"y0"=>[50.],"vx"=>[0.],"vy"=>[0.],"radius"=>[5.]),"tolerances"=>Dict("fx"=>NaN,"fy"=>NaN))
nothing # hide
```
## State and Control Names
```@example Bicycle
names = [:x,:y,:psi,:ux]
descriptions = ["X (m)","Y (m)","Yaw Angle (rad)","Longitudinal Velocity (m/s)"]
states!(n,names,descriptions=descriptions)
names = [:sa,:ax]
descriptions = ["Steering Angle (rad)","Longitudinal Acceleration (m/s^2)"]
controls!(n,names,descriptions=descriptions)
nothing # hide
```

## Differential Equations
```@example Bicycle
dx = KinematicBicycle_expr2(n)
dynamics!(n,dx)
nothing # hide
```

## Define and Configure the Problem:
```@example Bicycle
configure!(n;(:Nck=>[12,10,8]),(:finalTimeDV=>true));
nothing # hide
```

## Nonlinear Obstacle Avoidance Constraints
```@example Bicycle
sm = 2
x = n.r.ocp.x[:,1];y = n.r.ocp.x[:,2]; # pointers to JuMP variables
obs_con = @NLconstraint(n.ocp.mdl, [i=1:n.ocp.state.pts-1], 1 <= ((x[(i+1)]-c["obstacle"]["x0"][1])^2)/((c["obstacle"]["radius"][1]+sm)^2) + ((y[(i+1)]-c["obstacle"]["y0"][1])^2)/((c["obstacle"]["radius"][1]+sm)^2))
newConstraint!(n,obs_con,:obs_con)
nothing # hide
```

## Linear State Tolerances
```@example Bicycle
mXL = Any[false,false,false,false]
mXU = Any[false,false,false,-1];  # set to false if you don't want to taper that side
linearStateTolerances!(n;mXL=mXL,mXU=mXU)
nothing # hide
```

## Objective Function
```@example Bicycle
@NLobjective(n.ocp.mdl, Min, n.ocp.tf + (n.r.ocp.x[end,1]-c["goal"]["x"])^2 + (n.r.ocp.x[end,2]-c["goal"]["yVal"])^2);
nothing # hide
```

## Optimize
```@example Bicycle
optimize!(n);
nothing # hide
```

## Post Process
```@example Bicycle
plotSettings(;(:size=>(700,700)));
allPlots(n)
```
Notice the longitudinal velocity is pushed down to `29` m/s using the `linearStateTolerances!()` function.


The state limits can be turned off in the plots with `(:lims=>false)` and the obstacle plot handle can be passed to `statePlot()` in the 5th argument and by using `(:append=>true)`.


```@example Bicycle
plotSettings(;(:size=>(400,400)))
obs = obstaclePlot(n,c)
sp=statePlot(n,1,1,2,obs;(:append=>true),(:lims=>false))
xlims!(-45,55)
ylims!(0,110)
```
