using JuMP, Ipopt
Ni=2; Nck=[10,5];
τ = [zeros(Float64,Nck[int],) for int in 1:Ni];
τ[1] = [-1.0,-0.927484,-0.763842,-0.525646,-0.236234,0.0760592,0.380665,0.647767,0.851225,0.971175];
τ[2] = [-1.0,-0.72048,-0.167181,0.446314,0.885792];

mdl = Model(solver = IpoptSolver());

@variable(mdl, tf_var)

ts_JuMP = [Array(Any, Nck[int],) for int in 1:Ni];
for int in 1:Ni
  ts_JuMP[int] = tf_var/2*τ[int] + tf_var/2;
end

x = ts_JuMP[1];
N = length(x)
mat = repmat(x,1,N);

prod(mat,2)
