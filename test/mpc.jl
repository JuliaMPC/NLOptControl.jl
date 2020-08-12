using OrdinaryDiffEq, DiffEqBase

function MLplant(n,x0::Vector,t::Vector,U::Matrix,t0::Float64,tf::Float64)
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

```
integrationConfig=:lgrExplicit
interpolationOn=true
linearInterpolation=false
numInterpPts=250
AlltpolyPts=false
predictX0=false
maxSim=100
lastOptimal=false
tolOn=true
```
function MoonLander(integrationConfig;interpolationOn::Bool=true,linearInterpolation::Bool=false,numInterpPts::Int64=250,AlltpolyPts::Bool=false,predictX0::Bool=false,maxSim::Int64=100,lastOptimal::Bool=true,tolOn::Bool=false)
  n=define(numStates=2,numControls=1,X0=[10.,-2],XF=[0.,0.],CL=[0.],CU=[3.])
  n.s.ocp.numInterpPts = numInterpPts
  n.s.ocp.linearInterpolation = linearInterpolation
  n.s.ocp.interpolationOn = interpolationOn
  if tolOn
    X0_tol = [0.1,0.1]
    XF_tol = [0.1,0.1]
    defineTolerances!(n;X0_tol=X0_tol,XF_tol=XF_tol)
  end
  states!(n,[:h,:v];descriptions=["h(t)","v(t)"])
  controls!(n,[:T];descriptions=["T(t)"])
  dx=[:(v[j]),:(T[j]-1.625)]
  dynamics!(n,dx)
  configure!(n;(:integrationScheme=>integrationConfig),(:finalTimeDV=>true))
  obj=integrate!(n,:(T[j]))
  @NLobjective(n.ocp.mdl, Min, obj)
  initOpt!(n)
  defineMPC!(n;predictX0=predictX0,lastOptimal=lastOptimal,tex=0.1,printLevel=0,maxSim=maxSim)
  defineIP!(n,MLplant)
  simMPC!(n)

  if AlltpolyPts
    return (isequal(length(n.r.ocp.AlltpolyPts), length(n.r.ocp.dfs)) && n.f.mpc.goalReached)
  else
    return n.f.mpc.goalReached
  end
end


@testset "MoonLander MPC Tests" begin
  # Test #1
  @test MoonLander(:lgrExplicit)

  # Test #2
  #@test MoonLander(:lgrExplicit;lastOptimal=false,tolOn=true)

  # Test #3
  # a) :ps methods with several infeasibe solutions, but gets to the goal.
      # make sure that data lines up for plots even with infeasibe solutions
  # b) If using :ps methods and there is an infeasibe solution and the lastOptimal==true
      # if you try to use the interpolateLagrange!() function an error like this will occur:
        # WARNING: Not solved to optimality, status: Infeasible
        # ERROR: LoadError: BoundsError: attempt to access 39-element Array{Float64,1} at index [31:41]
      # To mitigate this error, if it was infeasibe and lastOptimal==true then a linearInterpolation is used.
  #@test MoonLander(:lgrExplicit;AlltpolyPts=true,numInterpPts=25,predictX0=true,maxSim=67)

  # Test #4
  #@test MoonLander(:lgrExplicit;linearInterpolation=true,predictX0=true,maxSim=72)

  # Test #5
  #@test MoonLander(:lgrExplicit;interpolationOn=false,predictX0=true)

  # Test #6
  #@test MoonLander(:bkwEuler;interpolationOn=false,predictX0=true)

  # Test #7
  #@test MoonLander(:trapezoidal;interpolationOn=false,predictX0=true)

  # Test #8
  # :tm methods with an infeasibe solution and a restoration error
 #@test MoonLander(:bkwEuler;predictX0=false)
end
