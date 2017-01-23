using NLOptControl
using JuMP
using Ipopt
using Parameters
using Plots
#using Gallium
#breakpoint_on_error()
#using VehicleModels
pyplot()

#TODO
# 1) why is it failing when we run it again?
# 4) try implicit method to see if that helps
# 6) try other linear solvers
# 9) allow for functions of state constraints of these so we can calculate them on the fly!
# 11) let the user define the objective function above
# 12) allow user to select from using the IMatrix or quadrature
# 13) make a bool to tell the user to restart julia if they change the model significantly -> or do this for them
# 14) make functionality to easily compare all intergration schemes

##################################
# Define NLOptControl problem
##################################

# Moon Lander Problem @ http://www.gpops2.com/Examples/MoonLander.html
const g = 1.62519; # m/s^2
# define dynamic constraint equations for constant final time
function stateEquations(x_int::Array{Any,2},u_int::Array{Any,2},st::Int64)
  if st==1
    return x_int[1:end-1,2]      # state eq# 1; v(t)
  elseif st==2
    return u_int[:,1] - g        # state eq# 2; u(t)-g
  end
end

#X0=[10.0,-2.0]; XF=[0.01,0.]
#XL=[-0.01,-Inf]; XU=[Inf,Inf];
#CL=[-Inf]; CU=[Inf];

X0=[10.0,-2.0]; XF=[0,0.]
XL=[-Inf,-Inf]; XU=[Inf,Inf];
CL=[-Inf]; CU=[Inf];

t0=0.0;tf=4.0;
 ps, nlp = initialize_NLP(numStates=2,
                          numControls=1,
                          Ni=1,Nck=[25],
                          stateEquations=stateEquations,
                          X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU;
                          (:finalTimeDV => true));

##################################
# Define JuMP problem
##################################
# initialize design problem
mdl = Model(solver = IpoptSolver(max_iter=3000)); #,warm_start_init_point = "yes"
#mdl = Model(solver=IpoptSolver(linear_solver = "mumps")) #linear_solver = "ma57"
d = JuMP.NLPEvaluator(mdl)
MathProgBase.initialize(d, [:Grad,:Hess, :Jac, :ExprGraph])

@unpack finalTimeDV = nlp
if finalTimeDV
  x,u,tf_var,ts_JuMP,ωₛ_JuMP,c = OCPdef(mdl,nlp,ps)
  obj = integrate(mdl,ps,u[:,1],ωₛ_JuMP;(:variable=>:control))
else
  x,u,c = OCPdef(mdl,nlp,ps)
  obj = integrate(mdl,ps,u[:,1];(:variable=>:control))
end

@NLobjective(mdl, Min, obj)
obj_val = solve(mdl)

#MathProgBase.constr_expr(d,i) #6
#MathProgBase.hesslag_structure(d)
#MathProgBase.jac_structure(d)
#MathProgBase.obj_expr(d)
##################################
# Post Processing
##################################
include(string(Pkg.dir("NLOptControl"),"/src/check_constraints.jl"))

if finalTimeDV; ts = ts_JuMP; else @unpack ts = ps; end
t_ctr= [idx for tempM in ts for idx = tempM[1:end-1]];
t_st = append!(t_ctr,ts[end][end]);

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
