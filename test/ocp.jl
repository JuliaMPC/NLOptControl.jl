BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig)) using expressions)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
  dynamics!(n,BrysonDenham_EXP)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,:(0.5*u1[j]^2));
  @NLobjective(n.ocp.mdl,Min,obj);
  optimize!(n);
  @show n.r.ocp.dfsOpt[:tSolve]
  @test isapprox(4/(9*L),n.r.ocp.objVal[1],atol=tol)
end

dx=[:(sin(x2[j])),:(u1[j])]
@testset "BeamProblem with (:integrationScheme=>$(integrationConfig)) using expressions and solver settings)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,XL=[-0.05,-1.0],XU=[-0.05,1.0]);
  dynamics!(n,dx)
  SS=[(:name=>:Ipopt),(:max_cpu_time=>5.0)]
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0),(:solverSettings=>SS));
  obj=integrate!(n,:( u1[j]^2 + 350*cos(x2[j]) ) )
  @NLobjective(n.ocp.mdl,Min,obj);
  optimize!(n);
  @show n.r.ocp.dfsOpt[:tSolve]
  @test isapprox(350,n.r.ocp.objVal[1],atol=tol)
end

dx=Array{Expr}(2);dx[1]=:(x[j]);dx[2]=:(u[j]-1.625);
@testset "MoonLander with (:integrationScheme=>$(integrationConfig)) using new names for states and controls; also naming them x and u to check for bugs)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.]);
  states!(n,[:h,:x];descriptions=["h(t)","x(t)"]);
  controls!(n,[:u];descriptions=["u(t)"]);
  dynamics!(n,dx)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>true));
  obj=integrate!(n,:(u[j]));
  @NLobjective(n.ocp.mdl, Min, obj);
  optimize!(n);
  @show n.r.ocp.dfsOpt[:tSolve]
  @test isapprox(8.9253,n.r.ocp.objVal[1],atol=tol)
end

############################
# Benchmarking Test
v0 = -2
h0 = 10
pts = 10000
tf_star = 2/3*v0 +4/3*(1/2*v0^2+3/2*h0)^0.5
ts = tf_star/2+v0/3

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

t_opt = linspace(0,tf_star,pts)
h_opt = zeros(pts); v_opt = zeros(pts); u_opt = zeros(pts);
for i in 1:pts
  h, v, u = optimal_solution(t_opt[i])
  h_opt[i] = h; v_opt[i] = v; u_opt[i] = u;
end

opt_runs = 2
col_pts = 10:2:12
Nck_vec = [[col_pts[i]] for i in 1:length(col_pts)]
opt_num = length(Nck_vec)

@testset "MoonLander test to make sure that a benchmark test will be able to be executed with (:integrationScheme=>$(integrationConfig)) " for integrationConfig in integrationConfigs

  # final average results
  t_opt_ave = NaN*zeros(opt_num)
  h_error_ave = NaN*zeros(opt_num)
  v_error_ave = NaN*zeros(opt_num)
  u_error_ave = NaN*zeros(opt_num)
  max_error_ave = NaN*zeros(opt_num)

  for num in 1:length(col_pts)
    n=define(numStates=2,numControls=1,X0=[h0,v0],XF=[0.,0.],XL=[-20,-20],XU=[20,20],CL=[0.],CU=[3.]);
    n.s.ocp.tfMax = 1000.0
    n.s.ocp.numInterpPts = pts  # must be the same as above to calculate error
    n.s.ocp.tfOptimal = tf_star # must be the same as above to calculate error
    states!(n,[:h,:v];descriptions=["h(t)","v(t)"]);
    controls!(n,[:T];descriptions=["T(t)"]);
    dx=[:(v[j]),:(T[j]-1.5)]
    dynamics!(n,dx)

    if (integrationConfig == :lgrImplicit || integrationConfig == :lgrExplicit)
      Nck = Nck_vec[num]
      configure!(n;(:integrationScheme=>integrationConfig),(:Nck=>Nck),(:finalTimeDV=>true));
    else
      configure!(n;(:integrationScheme=>integrationConfig),(:N=>col_pts[num]),(:finalTimeDV=>true));
    end
    x1=n.r.ocp.x[:,1];x2=n.r.ocp.x[:,2];
    obj=integrate!(n,:(T[j]));
    @NLobjective(n.ocp.mdl, Min, obj);
    setvalue(n.ocp.tf, 1.5)
    for i in 1:length(x1); setvalue(x1[i], 0.0); setvalue(x2[i], 0.0);  end
    # cache functions; inital optimization
    optimize!(n);

    # RESET temp variables for averaging results
    t_solve = NaN*zeros(opt_runs)
    h_error = NaN*zeros(opt_runs)
    v_error = NaN*zeros(opt_runs)
    u_error = NaN*zeros(opt_runs)
    max_error = NaN*zeros(opt_runs)

    for j in 1:opt_runs
      optimize!(n);
      if n.r.ocp.status == :Optimal
        t_solve[j] = n.r.ocp.tSolve
        h_error[j] = maximum(abs.(n.r.ocp.Xpts[:,1] - h_opt))
        v_error[j] = maximum(abs.(n.r.ocp.Xpts[:,2] - v_opt))
        u_error[j] = maximum(abs.(n.r.ocp.Upts[:,1] - u_opt))
        max_error[j] = maximum([h_error[j], v_error[j]])
      end
    end

    t_opt_ave[num] = mean(t_solve[.!isnan.(t_solve)],1)[1]
    h_error_ave[num] = mean(h_error[.!isnan.(h_error)],1)[1]
    v_error_ave[num] = mean(v_error[.!isnan.(v_error)],1)[1]
    u_error_ave[num] = mean(u_error[.!isnan.(u_error)],1)[1]
    max_error_ave[num] = mean(max_error[.!isnan.(max_error)],1)[1]
  end
  @show maximum(max_error_ave)
  @test isapprox(0, maximum(max_error_ave),atol=big_tol)
end

# Benchmarking Test
############################
