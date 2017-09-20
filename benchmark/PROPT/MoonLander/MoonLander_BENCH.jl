using NLOptControl,DataFrames,Interpolations,Plots

const v0 = -2
const h0 = 10
const ts = 1.4154
const tf_opt = 4.1641
const pts = 100

# actual optimal solution
function optimal_solution(t)
  if t < ts
    h = -3/4*t^2 + v0*t + h0
    v = -3/2*t + v0
    u = 0
  else
    h = 3/4*t^2 + (-3*ts + v0)*t + 3/2*ts^2 + h0
    v = 3/2*t + (-3*ts + v0)
    u = 3
  end
  return h, v, u
end
t_opt = linspace(0,tf_opt,pts)
h_opt = zeros(pts); v_opt = zeros(pts); u_opt = zeros(pts);

for i in 1:pts
  h, v, u = optimal_solution(t_opt[i])
  h_opt[i] = h; v_opt[i] = v; u_opt[i] = u;
end

scatter(t_opt,h_opt,label="h_opt")
scatter!(t_opt,v_opt,label="v_opt")
scatter!(t_opt,u_opt,label="u_opt")

# save the optimal solution for plotting and PROPT benchmark
dfs = DataFrame()
dfs[:t] = t_opt; dfs[:h] = h_opt; dfs[:v] = v_opt; dfs[:u] = u_opt;
writetable("optimal_solution.csv",dfs)

# define optimization in julia
opt_runs = 10
Nck_vec = [[60]]
opt_num = length(Nck_vec)

# final average results
t_opt_ave = zeros(opt_num)
h_error_ave = zeros(opt_num)
v_error_ave = zeros(opt_num)
u_error_ave = zeros(opt_num)
max_error_ave = zeros(opt_num)

# temp variables for averaging results
t_solve = zeros(opt_runs)
h_error = zeros(opt_runs)
v_error = zeros(opt_runs)
u_error = zeros(opt_runs)
max_error = zeros(opt_runs)
n=0;
for num in 1:opt_num
  Nck = Nck_vec[num]

  n=define(numStates=2,numControls=1,X0=[h0,v0],XF=[0.,0.],XL=[-20,-20],XU=[20,20],CL=[0.],CU=[3.]);
  n.s.tf_max=1000.0;
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"]);
  controls!(n,[:T];descriptions=["T(t)"]);
  dx=[:(v[j]),:(T[j]-1.5)]
  dynamics!(n,dx)
  configure!(n;(Nck=Nck),(:finalTimeDV=>true),(:solverSettings=>((:name=>:KNITRO),(:outlev=>0))));
  x1=n.r.x[:,1];x2=n.r.x[:,2];
  obj=integrate!(n,:(T[j]));
  @NLobjective(n.mdl, Min, obj);
  setvalue(n.tf, 1.5)
  for i in 1:Nck[1]; setvalue(x1[i], 0.0); setvalue(x2[i], 0.0);  end

  # cache functions; inital optimization
  optimize!(n);

  # temp variables for averaging results
  t_solve = zeros(opt_runs)
  h_error = zeros(opt_runs)
  v_error = zeros(opt_runs)
  u_error = zeros(opt_runs)
  max_error = zeros(opt_runs)

  for j in 1:opt_runs
    optimize!(n);

    # temporararily save results from this run
    t_solve[j] = n.r.t_solve

    # create splines TODO modify so :tm methods work here
    h_sp = interpolate((n.r.t_pts,),n.r.X_pts[:,1],Gridded(Linear()))
    v_sp = interpolate((n.r.t_pts,),n.r.X_pts[:,2],Gridded(Linear()))
    u_sp = interpolate((n.r.t_pts,),n.r.U_pts[:,1],Gridded(Linear()))

    h = zeros(pts); v = zeros(pts); u = zeros(pts);
    for i in 1:pts
      h[i] = h_sp[t_opt[i]]
      v[i] = v_sp[t_opt[i]]
      u[i] = u_sp[t_opt[i]]
    end

    h_error[j] = maximum(abs.(h - h_opt))
    v_error[j] = maximum(abs.(v - v_opt))
    u_error[j] = maximum(abs.(u - u_opt))
    max_error[j] = maximum([h_error[j], v_error[j], u_error[j]])
  end

  t_opt_ave[num] = mean(t_solve)
  h_error_ave[num] = mean(h_error)
  v_error_ave[num] = mean(v_error)
  u_error_ave[num] = mean(u_error)
  max_error_ave[num] = mean(max_error)
end

plot(h_error)
plot!(v_error)
plot!(u_error)

# optimize for a graph comparison






# save NLOptControl.jl results
#Nck
