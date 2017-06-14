

function BrysonDenham{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
  if n.s.integrationMethod==:tm; L=size(x)[1];else; L=size(x)[1]-1; end
  dx = Array(Any,L,n.numStates);
  dx[:,1] = @NLexpression(n.mdl,[j=1:L], x[j,2] );
  dx[:,2] = @NLexpression(n.mdl,[j=1:L], u[j,1] );
  return dx
end

@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig[1])) and (:integrationMethod=>$(integrationConfig[2]) using function)" for integrationConfig in integrationConfigs
  L=1/6;
  n=define!(BrysonDenham;numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN],CL=[NaN],CU=[NaN]);
  configure!(n;(:integrationMethod=>integrationConfig[1]),(:integrationScheme=>integrationConfig[2]),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,n.r.u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared));
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @test isapprox(4/(9*L),n.r.obj_val[1],atol=tol)
end

BrysonDenham_EXP=[:(x2[j]),:(u1[j])]

@testset "BrysonDenham with (:integrationScheme=>$(integrationConfig[1])) and (:integrationMethod=>$(integrationConfig[2]) using expressions)" for integrationConfig in integrationConfigs
  L=1/6;
  n=define!(BrysonDenham_EXP;numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[L,NaN],CL=[NaN],CU=[NaN]);
  configure!(n;(:integrationMethod=>integrationConfig[1]),(:integrationScheme=>integrationConfig[2]),(:finalTimeDV=>false),(:tf=>1.0));
  obj=integrate!(n,n.r.u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared));
  @NLobjective(n.mdl,Min,obj);
  optimize!(n);
  @test isapprox(4/(9*L),n.r.obj_val[1],atol=tol)
end
