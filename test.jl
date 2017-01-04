using NLOptControl
using JuMP
using Ipopt
using Parameters

# only want to call this once, even if repeatedly solving problem
ps, nlp = initialize_NLP(numStates=2,numControls=1,Ni=2,Nck=[3,4]);

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

# dynamic constraints
dζ = differentiate_state(ps,nlp)

# TODO check the length of the dynamic constraints-> should not be a dynamic constraint at the end
@NLconstraint(mdl, x_con[i=1:lengthStateVector-Ni], 0 == x[i+1] - x[i] )
@NLconstraint(mdl, v_con[i=1:lengthStateVector-Ni], 0 == v[i+1] - v[i] )

# defect constraints between intervals # TODO finish this
@constraint(mdl, x_defect[i=1:Ni-1], x[ ] == 0);
@constraint(mdl, v_defect[i=1:Ni-1], v[ ] == );

# calculate integral TODO need to integrate controls!
#ζ, approx_int_st = integrate_state(ps,nlp)  # for now, we can calculate the old way..

@NLobjective(mdl, Min, 0.5*u(t).^2 ) # TODO finish this

solve(mdl)
