using NLOptControl, JuMP, Parameters, PrettyPlots, Plots
pgfplots(); main_dir=pwd();

# initialize
n = NLOpt();

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
n = define(n,stateEquations=MoonLander,numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],XL=[-Inf,-Inf],XU=[Inf,Inf],CL=[0.],CU=[3.])

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
defineSolver(n,solver=:KNITRO)
XF_tol = [0.1, 0.1]; X0_tol = [0.1, 0.1]; defineTolerances(n;X0_tol=X0_tol,XF_tol=XF_tol);
names = [:h,:v]; descriptions = ["h(t)","v(t)"]; stateNames(n,names,descriptions);

# setup OCP
mdl = build(n);
n,r = OCPdef(mdl,n)
obj = integrate(mdl,n,r.u[:,1];C=1.0,(:variable=>:control),(:integrand=>:default))
@NLobjective(mdl, Min, obj);

# solve
optimize(mdl,n,r)

# post process
s=Settings(;format=:svg);
cd("results/")
allPlots(n,r,s,1)
cd(main_dir)
