using NLOptControl
using JuMP
using Ipopt
using Parameters
using Plots
#using Gallium
#breakpoint_on_error()
#using VehicleModels
pyplot()


##################################
# Define NLOptControl problem
##################################
n = NLOpt();
# Moon Lander Problem @ http://www.gpops2.com/Examples/MoonLander.html
const g = 1.62519; # m/s^2
function MoonLander{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2}) # dynamic constraint equations
  if n.integrationMethod==:tm
    L = size(x)[1];
  else
    L = size(x)[1]-1;
  end
  dx = Array(Any,L,n.numStates)
  dx[:,1] =  @NLexpression(mdl, [j=1:L], x[j,2] );
  dx[:,2] =  @NLexpression(mdl, [j=1:L], u[j,1] - g);
  return dx
end

n = define(n,stateEquations=MoonLander,numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],XL=[-Inf,-Inf],XU=[Inf,Inf],CL=[0.],CU=[3.])
#n = define(n,stateEquations=MoonLander,numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],XL=[-Inf,-Inf],XU=[Inf,Inf],CL=[-Inf],CU=[Inf])

#n = configure(n,Ni=2,Nck=[15,10];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 4.0))
#n = configure(n,N=10;(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => false),(:tf => 4.0))
#n = configure(n,N=10;(:integrationMethod => :tm),(:integrationScheme => :trapezoidal),(:finalTimeDV => false),(:tf => 4.0))
n = configure(n,Ni=3,Nck=[20,15,10];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV =>false),(:tf => 4.0))
#n = configure(n,N=30;(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => true))

##################################
# Define JuMP problem
##################################
# initialize design problem
mdl = Model(solver = IpoptSolver(max_iter=3000)); #,warm_start_init_point = "yes"

n,x,u,c=OCPdef(mdl,n)
obj = integrate(mdl,n,u[:,1];C=1.0,(:variable=>:control),(:integrand=>:default))

@NLobjective(mdl, Min, obj)
obj_val = solve(mdl)

##################################
# Post Processing
##################################
include(string(Pkg.dir("NLOptControl"),"/src/check_constraints.jl"))

if n.integrationMethod==:ps #TODO make this internal
  if n.finalTimeDV
    t_ctr= [idx for tempM in getvalue(n.ts) for idx = tempM[1:end-1]];
    t_st = append!(t_ctr,getvalue(n.ts[end][end]));
  else
    t_ctr= [idx for tempM in n.ts for idx = tempM[1:end-1]];
    t_st = append!(t_ctr,n.ts[end][end]);
  end
elseif n.integrationMethod==:tm
  if n.finalTimeDV
    t_ctr =  append!([0.0],cumsum(getvalue(n.dt)));
  else
    t_ctr =  append!([0.0],cumsum(n.dt));
  end
  t_st = t_ctr;
end

lw=8; lw2=3;
p1=plot(t_st,getvalue(x[:,1]), label = "x interp.",w=lw2)
scatter!(t_st,getvalue(x[:,1]), label = "x optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("h(t)")
xlabel!("time (s)")

p2=plot(t_st,getvalue(x[:,2]), label = "v interp.",w=lw2)
scatter!(t_st,getvalue(x[:,2]), label = "v optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("v(t)")
xlabel!("time (s)")

p3=plot(t_ctr,getvalue(u[:,1]), label = "u interp.",w=lw2)
scatter!(t_ctr,getvalue(u[:,1]), label = "u optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("u(t)")
xlabel!("time (s)")

plot(p1,p2,p3,layout=(1,3),background_color_subplot=RGB(0.2,0.2,0.2), background_color_legend=RGB(1,1,1))
plot!(foreground_color_grid=RGB(1,1,1))
