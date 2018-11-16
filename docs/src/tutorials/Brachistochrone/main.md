# Quick Ex#1: Brachistochrone
#### Solved by: John and Bernoulli, Newton and L'Hospital
### Given:
#### A particle sliding without friction along an unknown track in a gravitational field
#### Dynamic Constraints
$$\dot{x}_1(t)=x_3(t)sin(u_1(t))$$
$$\dot{x}_2(t)=-x_3(t)cos(u_1(t))$$
$$\dot{x}_3(t)=gcos(u_1(t))$$

#### Boundary Conditions
$${x}_1(0)=0 \qquad {x}_1(t_f)=2$$
$${x}_2(0)=0\qquad {x}_2(t_f)=-2$$
$${x}_3(0)=0\qquad {x}_3(t_f)=free$$
## Find:
#### The track that minimizes time
$$J=t_f$$
## Solution:

This problem can be found [here](http://www.gpops2.com/Examples/Brachistochrone.html).

## Packages that will be used
```@example Brachistochrone
using NLOptControl
nothing # hide
```

## Define the Problem
Next let's write down the boundary conditions into an array:
```@example Brachistochrone
X0=[0.0,0.0,0.0]
XF=[2.,-2.,NaN]
nothing # hide
```
#### Notice:
1. The numbers that where put into the expression are `Float64`; For now this is required!
2. The `NaN` is put into the boundary constraint for the third state; If any of the state bounds are free then pass a `NaN`

Now that we have the basic problem written down, we can call the `define()` function as:
```@example Brachistochrone
n=define(numStates=3,numControls=1,X0=X0,XF=XF);
nothing # hide
```
#### Basics
Variable | Description
:--- | :---
`n` | object that holds the entire optimal control problem
`de` | array of differential equation expressions
`numStates` | the number of states
`numControls` | the number of controls
`X0` | intial state constraint
`XF` | final state constraint

Also, not in this problem, but

Variable | Description
:--- | :---
`XL` | lower state bound
`XU` | upper state bound
`CL` | lower state bound
`CU` | upper control bound

The above bounds can be set in the same fashion as the initial and final state constraints. i.e. in an Array.

## State and Control Names (optional)
```@example Brachistochrone
states!(n,[:x,:y,:v],descriptions=["x(t)","y(t)","v(t)"]);
controls!(n,[:u],descriptions=["u(t)"]);
nothing # hide
```

## Differential Equations
Now we need to write all of the given information out. Let's start with the differential equation, that is written as an array of expressions:
```@example Brachistochrone
dx=[:(v[j]*sin(u[j])),:(-v[j]*cos(u[j])),:(9.81*cos(u[j]))]
dynamics!(n,dx)
nothing # hide
```

## Configure the Problem
Now that the basic optimal control problem has been defined, the next step is to `configure!()` with additional options.
```@example Brachistochrone
configure!(n;(:Nck=>[100]),(:finalTimeDV=>true));
nothing # hide
```
#### Settings:
Key | Description
:--- | :---
`:Nck` | array of that holds the number of points within each interval
`:finalTimeDV` | bool to indicate if time is a design variable

#### Notice:
1. Final time is a design variable; we are trying to minimize it
2. We defined this as a single interval problem with `100` points

## Objective Function
Finally, the objective function needs to be defined. For this, we use the `JuMP` macro `@NLOptControl()` directly as:
```@example Brachistochrone
@NLobjective(n.ocp.mdl,Min,n.ocp.tf)
nothing # hide
```
with,

Variable | Description
:--- | :---
`n.ocp.mdl` | object that holds them JuMP model
`Min` | for a minimization problem
`n.ocp.tf` | a reference to the final time

## Optimize
Now that the entire optimal control problem has been defined we can `optimize!()` it as:
```@example Brachistochrone
optimize!(n)
nothing # hide
```

## Post Process
**Make sure that you are not running the code in a folder where you have an important folder named `results`, because it will be deleted!**
Now that the problem has been optimized, we can quickly visualize the solution using ``allPlots()`` as:
```@example Brachistochrone
allPlots(n)
```
#### Optional plot settings
Many of the plot settings can be modified using the `plotSettings()` function. For instance;
```@example Brachistochrone
plotSettings(;(:size=>(700,700)));
```

`allPlots()` automatically plots the solution to all of the state and control variables. In this problem, we may be interested in comparing two states against one another which can be done using the `statePlot()` function as:
```@example Brachistochrone
statePlot(n,1,1,2)
```
For this case, there are four things that need to be passed to `statePlots()`:

Argument | Name | Description
:--- | :--- | :---
1|`n` |object that holds the entire optimal control problem
2|`idx` | reference to solution number used when we start solving `mpc` problems
3|`st1` | state number for `xaxis`
4|`st2` | state number for `yaxis`

#### Data Orginization
All of the states, control variables and time vectors are stored in an array of Dataframes called `n.r.dfs`
```@example Brachistochrone
n.r.ocp.dfs
```
It is an array because the problem is designed to be solved multiple times in a receding time horizon. The variables can be accessed like this:
```@example Brachistochrone
n.r.ocp.dfs[1][:x][1:4]
```

#### Optimization Data
```@example Brachistochrone
n.r.ocp.dfsOpt
```

The sailent optimization data is stored in the table above

Variable | Description
:--- | :---
`tSolve` | cpu time for optimization problem
`objVal` | objective function value
`iterNum` | a variable for a higher-level algorithm, often these problems are nested

One thing that may be noticed is the long time that it takes to solve the problem. This is typical for the first optimization, but after that even if the problem is modified the optimization time is greatly reduced.

#### For instance, let's re-run the optimization:
```@example Brachistochrone
optimize!(n);
n.r.ocp.dfsOpt[:tSolve]
```


## Costate visualization
For ``:ps`` methods the costates can also be calculates as

```@example Brachistochrone
X0=[0.0,0.0,0.0]
XF=[2.,-2.,NaN]
n=define(numStates=3,numControls=1,X0=X0,XF=XF)
states!(n,[:x,:y,:v],descriptions=["x(t)","y(t)","v(t)"])
controls!(n,[:u],descriptions=["u(t)"])
dx=[:(v[j]*sin(u[j])),:(-v[j]*cos(u[j])),:(9.81*cos(u[j]))]
dynamics!(n,dx)
n.s.ocp.evalCostates = true
configure!(n;(:Nck=>[100]),(:finalTimeDV=>true));
@NLobjective(n.ocp.mdl,Min,n.ocp.tf)
optimize!(n);
using PrettyPlots
allPlots(n)
```
Notice how the control jumps down for a bit, that is due to the equivalence of ``cos(n*2pi)`` for any integer `n`.

## Save results
While some results are save automatically, state, control, and costate (if applicable) data (about the collocation points and the Lagrange polynomial that runs through them) can be saved with the function `saveData()` as:
```@example Brachistochrone
saveData(n)
```
