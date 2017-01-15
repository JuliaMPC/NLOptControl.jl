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

# user defined dynamic constraint equations
function stateEquations(x_int::Array{Any,2},u_int::Array{Any,2},st::Int64)
  if st==1
    return x_int[1:end-1,2] # state eq# 1; v(t)
  elseif st==2
    return u_int[:,1]       # state eq# 2; u(t)
  end
end

L = 1/9;
X0=[0.;1]; XF=[0.;-1.] #TODO check to see what form this should be in ; or ,
                       #TODO allow for int inputs and just convert them to Float64
XL=[0.,-Inf]; XU=[L,Inf]; # TODO allow for functions of these so we can calculate them on the fly!
CL=[-Inf]; CU=[Inf];
ps, nlp = initialize_NLP(numStates=2,
                         numControls=1,
                         Ni=2,Nck=[5,5],
                         stateEquations=stateEquations,
                         X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU);

# given time interval--> not always given though; might be a design variable (especially for optimal control problems)
t0 = Float64(0); tf = Float64(1); @pack ps = t0, tf;

# give the time interval we can calculate these ps parameters
@unpack Nck, Ni, t0, tf, τ, ω = ps;
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω, ωₛ, ts;

##################################
# Define JuMP problem
##################################

# initialize design problem
mdl = Model(solver = IpoptSolver());

x,u=OCPdef(mdl,nlp,ps)

obj = integrate(mdl,ps,u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared))

@NLobjective(mdl, Min, obj)

obj_val = solve(mdl)

tol = 10e-3;
if abs(getvalue(obj) - 4/(9*L)) < tol
  print("\n Solution is correct to tolerance specs.!! \n \n")
else
  print(string("\n",
                "-------------------------------------", "\n",
                "The solution is not correct!!", "\n",
                "-------------------------------------", "\n",
                "The following values should be equal:", "\n",
                "4/(9*L)= ",4/(9*L),"\n",
                "getvalue(obj) = ",getvalue(obj),"\n"
                )
        )
end

####################################
## analytic soltion when 0<=L<=1/6
###################################
# analytic soltion for:    0 <= t <= 3L
u1(t) = -2/(3*L)*(1-t/(3*L));
v1(t) = (1 - t/(3*L))^2;
x1(t) = L*(1 - (1 - t/(3*L))^3 );

# analytic soltion for:   3L <= t <= 1-3L
u2(t) = 0;
v2(t) = 0;
x2(t) = L;

# analytic soltion for: 1-3L <= t <= 1
u3(t) = -2/(3*L)*(1 - (1-t)/(3*L) )
v3(t) = -(1 - (1-t)/(3*L) )^2
x3(t) =  L*(1 - (1 - (1-t)/(3*L) )^3 );

@unpack lengthStateVector = nlp
pts = 100;
t = linspace(t0,tf,pts)
x_analytic = zeros(Float64,pts,);
v_analytic = zeros(Float64,pts,);
u_analytic = zeros(Float64,pts,);

for i in 1:pts
  if L > 1/6
    warn("\n analytical solution only valid for L < 1/6!! \n")
  end

  if t[i] < 3*L
    u_analytic[i]=u1(t[i]);
    v_analytic[i]=v1(t[i]);
    x_analytic[i]=x1(t[i]);
  elseif ((3*L <= t[i]) && (t[i] <= (1-3*L)))
    u_analytic[i]=u2(t[i]);
    v_analytic[i]=v2(t[i]);
    x_analytic[i]=x2(t[i]);
  elseif (((1-3*L) <= t[i]) && (t[i] <= 1))
    u_analytic[i]=u3(t[i]);
    v_analytic[i]=v3(t[i]);
    x_analytic[i]=x3(t[i]);
  else
    error(" \n Not setup for outside of this range \n")
  end
end

# visualize solution
lw=8; lw2=3;
t_st = [idx for tempM in ts for idx = tempM];
t_ctr= [idx for tempM in ts for idx = tempM[1:end-1]];

p1=plot(t,x_analytic, label = "x analytic",w=lw)
plot!(t_st,getvalue(x[:,1]), label = "x interp.",w=lw2)
scatter!(t_st,getvalue(x[:,1]), label = "x optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("x(t)")
xlabel!("time (s)")

p2=plot(t,v_analytic, label = "v analytic",w=lw)
plot!(t_st,getvalue(x[:,2]), label = "v interp.",w=lw2)
scatter!(t_st,getvalue(x[:,2]), label = "v optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("v(t)")
xlabel!("time (s)")

p3=plot(t,u_analytic, label = "u analytic",w=lw)
plot!(t_ctr,getvalue(u[:,1]), label = "u interp.",w=lw2)
scatter!(t_ctr,getvalue(u[:,1]), label = "u optimal",marker = (:star8, 15, 0.9, :green))
ylabel!("u(t)")
xlabel!("time (s)")

plot(p1,p2,p3,layout=(1,3),background_color_subplot=RGB(0.2,0.2,0.2), background_color_legend=RGB(1,1,1))
plot!(foreground_color_grid=RGB(1,1,1))


#= COMMENTS old code and TODO stuff

#objective_fun(u1) = 0.5*integrate(u1,t0,tf);
#JuMP.register(:objective_fun,1,objective_fun,autodiff=true)

#=
#TODO add time as a optional design variable
@unpack numStatePoints, numControlPoints, XL, XU, CL, CU = nlp
@variable(mdl, XL[1] <= x[1:sum(numStatePoints)] <= XU[1])   # position
@variable(mdl, XL[2] <= v[1:sum(numStatePoints)] <= XU[2])   # velocity
@variable(mdl, CL[1] <= u[1:sum(numControlPoints)] <= CU[1]) # control #TODO make sure we can do LGR like this

#TODO automatically add constraints
# boundary constraints
@unpack X0, XF = nlp
@constraint(mdl, x0_con, x[1]   == X0[1]);
@constraint(mdl, xf_con, x[end] == XF[1]);
@constraint(mdl, v0_con, v[1]   == X0[2]);
@constraint(mdl, vf_con, v[end] == XF[2]);
=#


# TODO  allow for option to integrate constraints @unpack IMatrix
# TODO let the user define the objective function above
# TODO allow user to select from using the IMatrix or quadrature
=#
