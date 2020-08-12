# Bryson Denham


## Packages that will be used
using NLOptControl, PrettyPlots

n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[1/9,NaN]);

dx=[:(x2[j]),:(u1[j])]
dynamics!(n,dx)

configure!(n;(:integrationScheme=>:bkwEuler),(:finalTimeDV=>true));

obj=integrate!(n,:(0.5*u1[j]^2));
@NLobjective(n.ocp.mdl,Min,obj);

optimize!(n);
