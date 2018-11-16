# Quick Ex#2: Moon Lander

### Given:
#### A space-ship landing on the moon
#### Dynamic Constraints
$$\dot{x}_1(t)=x_2(t)$$
$$\dot{x}_2(t)=u(t)-g$$

#### Boundary Conditions
$${x}_1(0)=10 \qquad {x}_1(t_f)=0$$
$${x}_2(0)=-2 \qquad {x}_2(t_f)=0$$

#### Control Limits
$${u}_{1_{min}}=0$$
$${u}_{1_{max}}=3$$

### Find:
#### The track that minimizes time
$$J=\int_{0}^{tf} u(t) dt$$

## Solution:
In this problem, we put the bounds directly into `define()`. Also, now we have constant limits on the control variables and those can be added as shown below
This problem can be found [here](http://www.gpops2.com/Examples/MoonLander.html).

## Packages that will be used
```@example MoonLander
using NLOptControl
nothing # hide
```

## Define the Problem:
In this problem, we put the bounds directly into `define()`. Also, now we have constant limits on the control variables and those can be added as shown below
```@example MoonLander
n = define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
nothing # hide
```
## State and Control Names
The state and control variables are by default, $:x1,:x2,..$ and $:u1,:u2,..$, but they can be changed with the following commands:

```@example MoonLander
states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
controls!(n,[:T];descriptions=["T(t)"])
```
Next, now that the problem is configured, all of the state and control variables are stored in JuMP Arrays, `n.r.ocp.x[:,:]` and `n.r.u[:,:]`, respectively. For instance;
```@example MoonLander
typeof(n.r.ocp.x)
```
## Differential Equations
```@example MoonLander
dx=[:(v[j]),:(T[j]-1.625)]
dynamics!(n,dx)
nothing # hide
```

##  Configure the Problem:
```@example MoonLander
configure!(n;(:finalTimeDV=>true));
nothing # hide
```

## Integral Terms in the Cost Function
`integrate!()` is used to make terms that can be added to the cost function that need to be integrated. When calling this function an expression must be passed:

In this example the first control variable `T`needs to be integrated, it must be passed in an expression `:()` with the index `[j]`. To do this, `integrate!()` can be used as:
```@example MoonLander
obj = integrate!(n,:(T[j]))
# Now this term can be added as the objective function and the problem can be solved
@NLobjective(n.ocp.mdl, Min, obj);
nothing # hide
```

## Optimize
```@example MoonLander
optimize!(n)
nothing # hide
```

## Post Process
```@example MoonLander
allPlots(n)
```
## Other Dynamic Constraint Methods
Currently there are three different methods to ensure that the dyanamic constraints are satisfied and they are set when `configure!()` is called using the `:integrationScheme` key. They are listed below:

`:integrateScheme` | Description
:--- | :---
`:lgrExplicit`| default scheme; implementation derivative constraints in hp-pseudospecral method
`:lgrImplicit`| implementation of integral constraints in hp-pseudospecral method
`:bkwEuler` | approximate using backward euler method
`:trapezoidal` | approximate using trapezoidal method

The later two are time-marching methods and default number of points is `100`, but that can be changed by setting `N`. So, the above problem can be solved using one of the time-marching schemes as:
```@example MoonLander
n = define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
controls!(n,[:T];descriptions=["T(t)"])
dx=[:(v[j]),:(T[j]-1.625)]
dynamics!(n,dx)
configure!(n,N=200;(:integrationScheme=>:trapezoidal),(:finalTimeDV=>true))
obj = integrate!(n,:(T[j]))
@NLobjective(n.ocp.mdl, Min, obj)
optimize!(n)
allPlots(n)
```
## Constraints
Often when building a model and using it to solve an optimal control problem, their are issues associated with infeasibility. `NLOptControl` has functionality to help deal with these issues. For instance, the dual infeasibility values can stored and quickly viewed. They are stored in a DataFrame which can be referenced with `n.r.ocp.constraint.value` as:
```@example MoonLander
n.r.ocp.constraint.value
```
It is empty, because by default this data is not calculated and stored. This option can be turned on by modifying the settings for the problem:
```@example MoonLander
n.s.ocp.evalConstraints
```

```@example MoonLander
n.s.ocp.evalConstraints = true
optimize!(n)
n.r.ocp.constraint.value
```
```@example MoonLander
evalMaxDualInf(n)
```
The last function called, searches through all of the `dual infeasibilities` to find the largest value.
As, this problem is, it is feasible and optimal. But if there was an issue, often looking for high values in these DataFrame structures is the quickest way to figure out the constraints that are giving the solver trouble.

## Tolerances
If there was an example where the `dual infeasibility` value for one or more of the variables was very high, but the actual constraint is only being violated slightly (by some reasonable amount) then the tolerances on the initial and terminal states can be adjusted. This will also improve the solve time, so it is good practice to set these to reasonable values. For instance, in the `Moon Lander` example, we can set them as:
```@example MoonLander
n = define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
controls!(n,[:T];descriptions=["T(t)"])
dx = [:(v[j]),:(T[j]-1.625)]
dynamics!(n,dx)
XF_tol = [2.0,0.5]
X0_tol = [0.05,0.05]
defineTolerances!(n;X0_tol=X0_tol,XF_tol=XF_tol)
configure!(n,N=50;(:integrationScheme=>:bkwEuler),(:finalTimeDV=>true))
obj = integrate!(n,:(T[j]))
@NLobjective(n.ocp.mdl, Min, obj)
optimize!(n)
allPlots(n)
```
