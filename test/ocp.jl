BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig)) using expressions)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN]);
  dynamics!(n,BrysonDenham_EXP)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,:(0.5*u1[j]^2));
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(4/(9*L),n.r.obj_val[1],atol=tol)
end

dx=[:(sin(x2[j])),:(u1[j])]
@testset "BeamProblem with (:integrationScheme=>$(integrationConfig)) using expressions and solver settings)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,XL=[-0.05,-1.0],XU=[-0.05,1.0]);
  dynamics!(n,dx)
  SS=[(:name=>:Ipopt),(:max_cpu_time=>5.0)]
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0),(:solverSettings=>SS));
  obj=integrate!(n,:( u1[j]^2 + 350*cos(x2[j]) ) )
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(350,n.r.obj_val[1],atol=tol)
end

dx=Array{Expr}(2);
dx[1]=:(x[j])
dx[2]=:(u[j]-1.625)
@testset "MoonLander with (:integrationScheme=>$(integrationConfig)) using new names for states and controls; also naming them x and u to check for bugs)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.]);
  states!(n,[:h,:x];descriptions=["h(t)","x(t)"]);
  controls!(n,[:u];descriptions=["u(t)"]);
  dynamics!(n,dx)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>true));
  obj=integrate!(n,:(u[j]));
  @NLobjective(n.mdl, Min, obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(8.9253,n.r.obj_val[1],atol=tol)
end
