using NLOptControl
using JuMP
using Ipopt
using Parameters
using Plots

##################################
# Define NLOptControl problem
##################################

# user defined dynamic constraint equations
function stateEquations(nlp::NLP_data,ps::PS_data)
  @unpack ts, Ni, Nck, stateMatrix, controlMatrix = ps
  @unpack numStates = nlp
  dx = [zeros(Float64,length(ts[int])-1,numStates) for int in 1:Ni];
  for int in 1:Ni
    dx[int][:,1] = stateMatrix[int][:,2];   # state eq# 1; v(t) TODO either the stateMatrix needs to be updated or we need to put JuMP variables here
    dx[int][:,2] = controlMatrix[int][:,1]; # state eq# 2; u(t)
  end
  return dx
end
L = 1/9;
X0=[0.;1]; XF=[0.;-1.] #TODO check to see what form this should be in ; or ,
                       #TODO allow for int inputs and just convert them to Float64
XL=[0.,-Inf]; XU=[L,Inf]; # TODO allow for functions of these so we can calculate them on the fly!
CL=[-Inf]; CU=[Inf];
ps, nlp = initialize_NLP(numStates=2,
                         numControls=1,
                         Ni=1,Nck=[17],
                         stateEquations=stateEquations,
                         X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU);

# given time interval--> not always given though; might be a design variable (especially for optimal control problems)
t0 = Float64(0); tf = Float64(1);
@pack ps = t0, tf;

# give the time interval we can calculate these ps parameters
@unpack Nck, Ni, t0, tf, τ, ω = ps;
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω, ωₛ, ts;

##################################
# Define JuMP problem
##################################

# initialize design problem
mdl = Model(solver = IpoptSolver());

#TODO consider automatically generating JuMP variables based off user defined problem
#TODO add time as a optional design variable
# inequality constraints and design variable definition
@unpack numStatePoints, XL, XU, CL, CU = nlp
@variable(mdl, XL[1] <= x[1:sum(numStatePoints)] <= XU[1])   # position
@variable(mdl, XL[2] <= v[1:sum(numStatePoints)] <= XU[2])   # velocity
@variable(mdl, CL[1] <= u[1:sum(numStatePoints)] <= CU[1]) # control #TODO make sure we can do LGR like this

#TODO automatically add constraints
# boundary constraints
@unpack X0, XF = nlp
@constraint(mdl, x0_con, x[1]   == X0[1]);
@constraint(mdl, xf_con, x[end] == XF[1]);
@constraint(mdl, v0_con, v[1]   == X0[2]);
@constraint(mdl, vf_con, v[end] == XF[2]);

# state continuity constraints
@unpack Ni, Nck = ps;
Nck_st = [1;cumsum(Nck+1)];
for int in 2:Ni
  @constraint(mdl, x[Nck_st[int]] == x[Nck_st[int]+1])
  @constraint(mdl, v[Nck_st[int]] == v[Nck_st[int]+1])
  @constraint(mdl, u[Nck_st[int]] == u[Nck_st[int]+1])
end

# calculate LGR matrices - > IMatrix and DMatrix
LGR_matrices(ps,nlp); # TODO if the final time is changing, DMatrix will change!--> needs to be recalcualted during optimization
@unpack DMatrix, t0, tf = ps; # TODO look into regestering the DMatrix

# state equation constraints
Nck_st  = [0;cumsum(Nck+1)];
#Nck_ctr = [0;cumsum(Nck)];
for int in 1:Ni
  # states
  x_int = x[Nck_st[int]+1:Nck_st[int+1]];
  v_int = v[Nck_st[int]+1:Nck_st[int+1]];

  # controls
  u_int = u[Nck_st[int]+1:Nck_st[int+1]];

  # dynamics
  dx1 = DMatrix[int]*x_int - v_int[1:end-1];# TODO automatically insert user defined stateEquations()
  dx2 = DMatrix[int]*v_int - u_int[1:end-1];

  # add dynamics constraints
  for j in 1:Nck[int]
    @constraint(mdl, 0 == dx1[j])
    @constraint(mdl, 0 == dx2[j])
  end
end

# TODO  allow for option to integrate constraints @unpack IMatrix

# TODO let the user define the objective function above
# TODO allow user to select from using the IMatrix or quadrature

#u_int = u[Nck_ctr[int]+1:Nck_ctr[int+1]];
#TODO find a better way to do this!
@unpack ωₛ = ps
@NLexpression(mdl, temp[int=1:Ni], 0.5*sum{ωₛ[int][j]*u[Nck_ctr[int]+1:Nck_ctr[int+1]][j]*u[Nck_ctr[int]+1:Nck_ctr[int+1]][j],j=1:Nck[int]});
@NLexpression(mdl, obj, sum{temp[int], int=1:Ni});
# TODO check obj to make sure that it is not an array

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

  #TODO why is control going to zero between intervals?
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

t = linspace(t0,tf,lengthStateVector)  #TODO interpolate using Lagrange Polynomial?
u_analytic = zeros(Float64,lengthStateVector,);
v_analytic = u_analytic;
x_analytic = u_analytic;
for i in 1:lengthStateVector
  if L > 1/6
    warn("\n analytical solution only valid for L < 1/6!! \n")
  end

  if t[i] < 3*L
    u_analytic[i]=u1(t[i]);
    v_analytic[i]=v1(t[i]);
    x_analytic[i]=x1(t[i]);
    print(1)
  elseif ((3*L <= t[i]) && (t[i] <= (1-3*L)))
    u_analytic[i]=u2(t[i]);
    v_analytic[i]=v2(t[i]);
    x_analytic[i]=x2(t[i]);
    print(2)
  elseif (((1-3*L) <= t[i]) && (t[i] <= 1))
    u_analytic[i]=u3(t[i]);
    v_analytic[i]=v3(t[i]);
    x_analytic[i]=x3(t[i]);
    print(3)
  else
    error(" \n not setup for outside of this range \n")
  end
end

# visualize solution TODO plot analytic solution
lw=4; lw2=3;
tdata = [idx for tempM in ts for idx = tempM];

p1=plot(t,x_analytic, label = "x analytic",w=lw)
plot!(tdata,getvalue(x), label = "x interp.",w=lw2)
scatter!(tdata,getvalue(x), label = "x optimal")
ylabel!("x(t)")
xlabel!("time (s)")

p2=plot(t,v_analytic, label = "v analytic",w=lw)
plot!(tdata,getvalue(v), label = "v interp.",w=lw2)
scatter!(tdata,getvalue(v), label = "v optimal")
ylabel!("v(t)")
xlabel!("time (s)")

p3=plot(t,u_analytic, label = "u analytic",w=lw)
plot!(tdata,getvalue(u), label = "u interp.",w=lw2)
scatter!(tdata,getvalue(u), label = "u optimal")
ylabel!("u(t)")
xlabel!("time (s)")

plot(p1,p2,p3,layout=(1,3))
