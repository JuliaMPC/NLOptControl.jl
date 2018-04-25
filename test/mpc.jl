using OrdinaryDiffEq, DiffEqBase


function IPplant(n,x0::Vector,t::Vector,U::Matrix,t0::Float64,tf::Float64)
  spU = linearSpline(t,U[:,1])   # create spline
  f = (dx,x,p,t) -> begin        # diff eqs.
    dx[1] = x[2]
    dx[2] = spU[t] - 1.625
  end
  tspan = (t0,tf)
  prob = ODEProblem(f,x0,tspan)
  sol = DiffEqBase.solve(prob,Tsit5())
  U = [spU]
  return sol, U
end

@testset "MPC Tests" begin
  # MoonLander Test #1
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
  controls!(n,[:T];descriptions=["T(t)"])
  dx=[:(v[j]),:(T[j]-1.625)]
  dynamics!(n,dx)
  configure!(n;(:finalTimeDV=>true))
  obj=integrate!(n,:(T[j]))
  @NLobjective(n.ocp.mdl, Min, obj)
  initOpt!(n)
  defineMPC!(n;predictX0=false,tex=0.1,printLevel=0)
  defineIP!(n,IPplant)
  simMPC!(n)
  @test n.f.mpc.goalReached

  # MoonLander test #2
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
  X0_tol = [0.1,0.1]
  XF_tol = [0.1,0.1]
  defineTolerances!(n;X0_tol=X0_tol,XF_tol=XF_tol)
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
  controls!(n,[:T];descriptions=["T(t)"])
  dx=[:(v[j]),:(T[j]-1.625)]
  dynamics!(n,dx)
  configure!(n;(:finalTimeDV=>true))
  obj=integrate!(n,:(T[j]))
  @NLobjective(n.ocp.mdl, Min, obj)
  initOpt!(n)
  defineMPC!(n;predictX0=false,tex=0.1,lastOptimal=false,printLevel=0)
  defineIP!(n,IPplant)
  simMPC!(n)
  @test n.f.mpc.goalReached

  # MoonLander Test #3  (:tm methods with an infeasibe solution and a restoration error)
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
  controls!(n,[:T];descriptions=["T(t)"])
  dx=[:(v[j]),:(T[j]-1.625)]
  dynamics!(n,dx)
  configure!(n;(:integrationScheme=>:bkwEuler),(:finalTimeDV=>true))
  obj=integrate!(n,:(T[j]))
  @NLobjective(n.ocp.mdl, Min, obj)
  initOpt!(n)
  defineMPC!(n;predictX0=false,tex=0.1,printLevel=0)
  defineIP!(n,IPplant)
  simMPC!(n)
  @test n.f.mpc.goalReached

end
