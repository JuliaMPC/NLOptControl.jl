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

de=[:(sin(x2[j])),:(u1[j])]
@testset "BeamProblem with (:integrationScheme=>$(integrationConfig)) using expressions)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,XL=[-0.05,-1.0],XU=[-0.05,1.0]);
  dynamics!(n,de)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,:( u1[j]^2 + 350*cos(x2[j]) ) )
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(350,n.r.obj_val[1],atol=tol)
end

de=[:(v[j]),:(T[j]-1.625)]
@testset "MoonLander with (:integrationScheme=>$(integrationConfig)) using new names for states and controls)" for integrationConfig in integrationConfigs
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.]);
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"]);
  controls!(n,[:T];descriptions=["T(t)"]);
  dynamics!(n,de)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>true));
  obj=integrate!(n,:(T[j]));
  @NLobjective(n.mdl, Min, obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(8.9253,n.r.obj_val[1],atol=tol)
end
