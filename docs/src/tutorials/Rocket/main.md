# Rocket

This problem can be found [here](https://github.com/JuliaOpt/juliaopt-notebooks/blob/master/notebooks/JuMP-Rocket.ipynb).

## Packages that will be used
```@example Rocket
using NLOptControl
nothing # hide
```

## Constants
```@example Rocket
# Note that all parameters in the model have been normalized
# to be dimensionless. See the COPS3 paper for more info.
h_0 = 1    # Initial height
v_0 = 0    # Initial velocity
m_0 = 1    # Initial mass
g_0 = 1    # Gravity at the surface

# Parameters
T_c = 3.5  # Used for thrust
h_c = 500  # Used for drag
v_c = 620  # Used for drag
m_c = 0.6  # Fraction of initial mass left at end

# Derived parameters
c     = 0.5*sqrt(g_0*h_0)  # Thrust-to-fuel mass
m_f   = m_c*m_0            # Final mass
D_c   = 0.5*v_c*m_0/g_0    # Drag scaling
T_max = T_c*g_0*m_0        # Maximum thrust
nothing # hide
```

## Define the Problem:
```@example Rocket
n=define(numStates=3,numControls=1,X0=[h_0,v_0,m_0],XF=[NaN,NaN,m_f],XL=[h_0,v_0,m_f],XU=[NaN,NaN,m_0],CL=[0.0],CU=[T_max]);
nothing # hide
```

## State and Control Names
```@example Rocket
states!(n,[:h,:v,:m],descriptions=["height (t)","velocity (t)","mass (t)"]);
controls!(n,[:T],descriptions=["thrust (t)"]);
```

## Differential Equations
```@example Rocket
Drag=:($D_c*v[j]^2*exp(-$h_c*(h[j]-$h_0)/$h_0));
Grav=:($g_0*($h_0/h[j])^2);
dx=Array{Expr}(3,);
dx[1]=:(v[j]);
dx[2]=:((T[j]-$Drag)/m[j]-$Grav)
dx[3]=:(-T[j]/$c);
```
Then add the differential equations to the model:
```@example Rocket
dynamics!(n,dx)
```

## Configure the Problem:
```@example Rocket
configure!(n;(:finalTimeDV=>true));
nothing # hide
```

## Objective Function
```@example Rocket
@NLobjective(n.ocp.mdl,Max,n.r.ocp.x[end,1]);
nothing # hide
```

## Optimize
```@example Rocket
optimize!(n);
nothing # hide
```

## Post Process
```@example Rocket
allPlots(n)
```
