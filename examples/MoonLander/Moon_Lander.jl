using NLOptControl, JuMP, Parameters
main_dir=pwd();s=Settings();

n = NLOpt(); # initialize

# Moon Lander Problem @ http://www.gpops2.com/Examples/MoonLander.html
const g = 1.62519; # m/s^2
function MoonLander{T<:Any}(mdl::JuMP.Model,n::NLOpt,r::Result,x::Array{T,2},u::Array{T,2}) # dynamic constraint equations
  if n.integrationMethod==:tm; L=size(x)[1]; else L=size(x)[1]-1; end
  dx = Array(Any,L,n.numStates)
  dx[:,1] =  @NLexpression(mdl, [j=1:L], x[j,2] );
  dx[:,2] =  @NLexpression(mdl, [j=1:L], u[j,1] - g);
  return dx
end

# define
#n = define(n,stateEquations=MoonLander,numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],XL=[NaN,NaN],XU=[NaN,NaN],CL=[0.],CU=[3.])
n = define(n,stateEquations=MoonLander,numStates=2,numControls=1,X0=[10.,-2],XF=[0.,NaN],XL=[NaN,-3],XU=[NaN,-1],CL=[0.],CU=[3.])

# build
# no time
#n = configure(n,Ni=2,Nck=[15,10];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 4.0))
#n = configure(n,N=10;(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => false),(:tf => 4.0))
#n = configure(n,N=10;(:integrationMethod => :tm),(:integrationScheme => :trapezoidal),(:finalTimeDV => false),(:tf => 4.0))

# with time
n = configure(n,Ni=4,Nck=[5,5,4,6];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV =>true))
#n = configure(n,N=30;(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => true))
#n = configure(n,N=30;(:integrationMethod => :tm),(:integrationScheme => :trapezoidal),(:finalTimeDV => true))

# addtional information
mXL=Any[false,0.1];mXU=Any[false,-.1];  # set to false if you don't want to taper that side
linearStateTolerances(n;mXL=mXL,mXU=mXU);
#defineSolver(n,solver=:KNITRO)
#XF_tol = [0.001, 0.001]; X0_tol = [0.001, 0.001]; defineTolerances(n;X0_tol=X0_tol,XF_tol=XF_tol);
names = [:h,:v]; descriptions = ["h(t)","v(t)"]; stateNames(n,names,descriptions);

# setup OCP
mdl = build(n);
n,r = OCPdef(mdl,n,s)
obj = integrate(mdl,n,r.u[:,1];C=1.0,(:variable=>:control),(:integrand=>:default))
@NLobjective(mdl, Min, obj);

# solve
optimize(mdl,n,r,s)

# post process
using PrettyPlots, Plots
gr();
s=Settings(;format=:png);
allPlots(n,r,s,1)
