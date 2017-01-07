using NLOptControl
using JuMP
using Ipopt
using Parameters

# only want to call this once, even if repeatedly solving problem
ps, nlp = initialize_NLP(numStates=2,numControls=1,Ni=1,Nck=[4]);

# given time interval--> not always given though; might be a design variable (especially for optimal control problems)
t0 = Float64(0); tf = Float64(1);
@pack ps = t0, tf;

# give the time interval we can calculate these ps parameters
@unpack Nck, Ni, t0, tf, τ, ω = ps;
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω, ωₛ, ts;

# now the LGR matrices can be calculated
LGR_matrices(ps,nlp) # will have to do this inside the optimization if the final time is a design variable

# needed to define lengths of dvs
@unpack lengthStateVector, lengthControlVector = nlp

# initialize design problem
mdl = Model(solver = IpoptSolver());

# design variables
@variable(mdl, x[1:lengthStateVector] <= 1/9) # position
@variable(mdl, v[1:lengthStateVector])        # velocity
@variable(mdl, u[1:lengthControlVector])      # control

# boundary constraints
@constraint(mdl, x0_con, x[1] == 0);
@constraint(mdl, xf_con, x[end] == 0);
@constraint(mdl, v0_con, v[1] == 1);
@constraint(mdl, vf_con, v[end] == -1);

# define dynamic constraint equations
dx(t) = v(t);
dv(t) = u(t);

dx(t) = v[t];
dv(t) = u[t];

function ODE_constraint(pa::PS_data,nlp::NLP_data)
@unpack stateMatrix, controlMatrix



return
end
if p.scale == 1
    x = scaling(x,p.xL,p.xU,2); # unscale design variables
end

# initialize constraints
c = []; ceq = []; contin = [];
# extract from optimization vector
X = reshape(x(1:p.ns*p.nt),p.nt,[]); # states
U = reshape(x(p.ns*p.nt+1:end),p.nt,[]); # controls

#--- equality constraints
# determine defect constraints for each segment

#temp
#    p.nt = sum(p.Narray); # total number of nodes
#    p.Narray = Nc
#    p.cumN = [0,cumsum(p.Narray)]; # segment ending indices
#temp

for i = 1:length(p.Narray)
    Xi = X(p.cumN(i)+1:p.cumN(i+1),:); # states for segment i
    Ui = U(p.cumN(i)+1:p.cumN(i+1),:); # controls for segment i
    D = p.D{i}; # differentiation matrix

    for i = 1:N+1
        f(i,:) = eval([p.deriv,'(tau(i),X(i,:),U(i,:),t0,tf,p)']);
    end
    F = (tf-t0)/2*f; % scale

    F = Fmatrix(p.tau{i},Xi,Ui,p.t0(i),p.tf(i),p); # differential matrix
    defect{i,1} = D*Xi - F; # defect constraints
end
ceq = [ceq;reshape(cell2mat(defect),[],1)];

# state and control continuity constraints between segments
for i = 1:length(p.Narray)-1
    contin{i,1} = [X(p.cumN(i+1),:) - X(p.cumN(i+1)+1,:),...
        U(p.cumN(i+1),:) - U(p.cumN(i+1)+1,:)];
end
ceq = [ceq;reshape(cell2mat(contin),[],1)];

# boundary constraints
boundary = eval([p.boundary,'(X,U,p)']);
ceq = [ceq;reshape(boundary,[],1)];

#--- inequality constraints
# path constraints
path = eval([p.path,'(X,U,p)']);
c = [c;reshape(path,[],1)];





# dynamic constraints
dζ = differentiate_state(ps,nlp)

@unpack DMatrix
# TODO check the length of the dynamic constraints-> should not be a dynamic constraint at the end
# TODO are these linear or nonlinear constraints?
# TODO brush up on sum fucntion in JuMP
# TODO can we do tripple for loops in JuMP constraints?
@constraint(mdl, x_con[j=1:Nck[int]+1], 0 == sum( DMatrix[int][i][j]*x[j] - (tf - t0)/2*v[i]) )
@constraint(mdl, v_con[i=1:lengthStateVector-Ni], 0 == v[i+1] - v[1] )

# TODO  allow for Integration of constraints
if false
  @unpack IMatrix
  @NLconstraint(mdl, x_con[i=1:lengthStateVector-Ni], 0 == x[i+1] - x[1] - (tf - t0)/2* I??)
  @NLconstraint(mdl, v_con[i=1:lengthStateVector-Ni], 0 == v[i+1] - v[1] )
end
# defect constraints between intervals # TODO finish this
@constraint(mdl, x_defect[i=1:Ni-1], x[ ] == 0);
@constraint(mdl, v_defect[i=1:Ni-1], v[ ] == );

# calculate integral TODO need to integrate controls!
#ζ, approx_int_st = integrate_state(ps,nlp)  # for now, we can calculate the old way..

@NLobjective(mdl, Min, 0.5*u(t).^2 ) # TODO finish this

solve(mdl)
