

function BrysonDenham{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
  if n.s.integrationMethod==:tm; L=size(x)[1];else; L=size(x)[1]-1; end
  dx = Array(Any,L,n.numStates);
  dx[:,1] = @NLexpression(n.mdl,[j=1:L], x[j,2] );
  dx[:,2] = @NLexpression(n.mdl,[j=1:L], u[j,1] );
  return dx
end

L=1/6;
@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig)) using function)" for integrationConfig in integrationConfigs
  n=define!(BrysonDenham;numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN],CL=[NaN],CU=[NaN]);
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,n.r.u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared));
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(4/(9*L),n.r.obj_val[1],atol=tol)
end

BrysonDenham_EXP=[:(x2[j]),:(u1[j])]; L=1/6;
@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig)) using expressions)" for integrationConfig in integrationConfigs
  n=define!(BrysonDenham_EXP;numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN],CL=[NaN],CU=[NaN]);
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,n.r.u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared));
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(4/(9*L),n.r.obj_val[1],atol=tol)
end


de=[:(sin(x2[j])),:(u1[j])]
@testset "BeamProblem with (:integrationScheme=>$(integrationConfig)) using expressions)" for integrationConfig in integrationConfigs
  n=define!(de;numStates=2,numControls=1,X0=[NaN,NaN],XF=[NaN,NaN],XL=[-0.05,-1.0],XU=[-0.05,1.0],CL=[NaN],CU=[NaN]);
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>false),(:tf=>1.0));
  obj1=integrate!(n,n.r.u[:,1];(:variable=>:control),(:integrand=>:squared));
  obj2=integrate!(n,n.r.x[:,2];C=350.,(:variable=>:state),(:integrand=>:cos));
  @NLobjective(n.mdl,Min,obj1+obj2);
  optimize!(n);
  @show n.r.dfs_opt[1][:t_solve][1]
  @test isapprox(350,n.r.obj_val[1],atol=tol)
end
