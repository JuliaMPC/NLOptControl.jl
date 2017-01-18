using NLOptControl
using JuMP
using Ipopt
using Parameters
using Plots
#using VehicleModels
pyplot()

##################################
# Define NLOptControl problem
##################################

# Moon Lander Problem
# http://www.gpops2.com/Examples/MoonLander.html
const g = 1.62519; # m/s^2
# define dynamic constraint equations
function stateEquations(mdl::JuMP.Model,x_int::Array{Any,2},u_int::Array{Any,2},st::Int64)
  if st==1
    return x_int[1:end-1,2]      # state eq# 1; v(t)
  elseif st==2
    return u_int[:,1] - g        # state eq# 2; u(t)-g
  end
end

X0=[10.0,-2.0]; XF=[0.,0.]    # TODO allow for int inputs and just convert them to Float64
XL=[-Inf,-Inf]; XU=[Inf,Inf]; # TODO allow for functions of these so we can calculate them on the fly!
CL=[-Inf]; CU=[Inf];
#t0 = Float64(0); tf = Float64(4);
ps, nlp = initialize_NLP(numStates=2,
                         numControls=1,
                         Ni=2,Nck=[10,5],
                         stateEquations=stateEquations,
                         X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU;
                         (:finalTimeDV => true));

#=
t0 = Float64(0); tf = Float64(4); @pack ps = t0, tf;
# give the time interval we can calculate these ps parameters
@unpack Nck, Ni, t0, tf, τ, ω = ps;
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω, ωₛ, ts;
=#

##################################
# Define JuMP problem
##################################

# initialize design problem
#mdl = Model(solver = IpoptSolver());
mdl = Model(solver=IpoptSolver(linear_solver = "mumps")) #linear_solver = "ma57"

d = JuMP.NLPEvaluator(mdl)
MathProgBase.initialize(d, [:Grad,:Hess, :Jac, :ExprGraph])
#MathProgBase.hesslag_structure(d)
#MathProgBase.jac_structure(d)
 #MathProgBase.obj_expr(d)
 for i in 1:20
  print(i)
 MathProgBase.constr_expr(d,7) #6
 end
x,u = OCPdef(mdl,nlp,ps)

obj = integrate(mdl,ps,u[:,1];(:variable=>:control))

@NLobjective(mdl, Min, obj) # for now might actually be able to remove this! tf_var

obj_val = solve(mdl)

##################################
# Post Processing
##################################
@unpack ts = ps
lw=8; lw2=3;
t_ctr= [idx for tempM in ts for idx = tempM[1:end-1]];
t_st = append!(t_ctr,ts[end][end]);

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
