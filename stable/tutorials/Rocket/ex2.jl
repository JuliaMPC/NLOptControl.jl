# Rocket
using NLOptControl, PrettyPlots
n = define(numStates=3,numControls=1,X0=[0.0,0.0,1.0],XF=[10.0,0.0,NaN],XL=[NaN,-0.1,0.0],XU=[NaN,1.7,1.0],CL=[-1.1],CU=[1.1]);
n.s.ocp.tfMax = 15.; # tFmin = 5.
states!(n,[:h,:v,:m],descriptions=["height (t)","velocity (t)","mass (t)"]);
controls!(n,[:T],descriptions=["thrust (t)"]);
dx=Array{Expr}(3,);
dx[1] = :(v[j])
dx[2] = :((T[j]- 0.2*v[j]^2)/m[j])
dx[3] = :(-0.01*T[j]^2)
dynamics!(n,dx)
configure!(n;(:solverSettings=>(:name=>:KNITRO)),(:finalTimeDV=>true))
@NLobjective(n.ocp.mdl,Min, n.ocp.tf)
optimize!(n)
allPlots(n)
