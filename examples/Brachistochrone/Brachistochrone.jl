using NLOptControl, JuMP, Parameters, Plots, KNITRO
gr(); main_dir=pwd();
n = NLOpt(); s=Settings(); # initialize

# Brachistichrone Problem @ http://gpops2.com/Examples/Brachistochrone.html
                        # @ http://apmonitor.com/wiki/index.php/Apps/BrachistochroneProblem
const g = 9.81; # m/s^2
function Brachistichrone{T<:Any}(mdl::JuMP.Model,n::NLOpt,r::Result,x::Array{T,2},u::Array{T,2}) # dynamic constraint equations
  if n.integrationMethod==:tm; L=size(x)[1]; else L=size(x)[1]-1; end
  dx = Array(Any,L,n.numStates)
  dx[:,1] =  @NLexpression(mdl, [j=1:L], x[j,3]*sin(u[j,1]) );
  dx[:,2] =  @NLexpression(mdl, [j=1:L], x[j,3]*cos(u[j,1]) );
  dx[:,3] =  @NLexpression(mdl, [j=1:L], g*cos(u[j,1]) );
  return dx
end

# define
n = define(n,stateEquations=Brachistichrone,numStates=3,numControls=1,X0=[0.0,0.0,0.0],XF=[2.,2.,NaN],XL=[-Inf,-Inf,-Inf],XU=[Inf,Inf,Inf],CL=[-Inf],CU=[Inf])

n = configure(n,Ni=4,Nck=[5,5,4,6];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV =>true))
#n = configure(n,N=30;(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => true))
#n = configure(n,N=30;(:integrationMethod => :tm),(:integrationScheme => :trapezoidal),(:finalTimeDV => true))

# addtional information
names = [:x,:y,:v]; descriptions = ["x(t)","y(t)","v(t)"]; stateNames(n,names,descriptions);
names = [:u]; descriptions = ["u(t)"]; controlNames(n,names,descriptions);

# setup OCP
defineSolver(n,solver=:KNITRO)  # TODO figure out why this does not run using IPOPT!!
mdl = build(n);
n,r = OCPdef(mdl,n)
@NLobjective(mdl, Min, n.tf);

optimize(mdl,n,r,s) # solve

# post process
s=Settings(;format=:png);
allPlots(n,r,s,1)
