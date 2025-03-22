"""
"""
function intervals(n::NLOpt,int,x,u)

  if typeof(x[1,1]) == JuMP.Variable

    # States
    x_int = Array{JuMP.Variable}(length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.state.num)
    for st in 1:n.ocp.state.num # +1 adds the DV in the next interval
      x_int[:,st] = x[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,st]
    end

    # Controls
    if int!=n.ocp.Ni
      u_int = Matrix{JuMP.Variable}(undef, length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.control.num)
    else                    # -1 -> removing control in last mesh interval
      u_int = Matrix{JuMP.Variable}(undef, length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]-1),n.ocp.control.num)
    end
    for ctr in 1:n.ocp.control.num
      if int!=n.ocp.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,ctr]
      else
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1],ctr]
      end
    end

  else
    # states
    x_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.state.num);
    for st in 1:n.ocp.state.num # +1 adds the DV in the next interval
      x_int[:,st] = x[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,st]
    end

    # controls
    if int!=n.ocp.Ni
      u_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]),n.ocp.control.num)
    else                    # -1 -> removing control in last mesh interval
      u_int = Matrix{Any}(undef,length(n.ocp.Nck_full[int]+1:n.ocp.Nck_full[int+1]-1),n.ocp.control.num)
    end
    for ctr in 1:n.ocp.control.num
      if int!=n.ocp.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1,ctr]
      else
        u_int[:,ctr] = u[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1],ctr]
      end
    end
  end

  return x_int,u_int
end

"""
"""
function interpolateLagrange!(n::NLOpt{T}; numPts::Int=250, tfOptimal::Any=false) where { T <: Number }
  # TODO throw an error if tfOptimal does not make sense given current solution
  if isa(tfOptimal,Bool)
    if n.s.ocp.finalTimeDV
      tf = getvalue(n.ocp.tf)
    else n.s.ocp.finalTimeDV
      tf = n.ocp.tf
    end
  else  # if there is a known optimal final time, then it is useful to evaluate the Lagrange polynomial at particular points to determine the error in the solution
    tf = tfOptimal
  end

  if isnan(tf)
      @warn "tf is a NaN cannot use it in interpolateLagrange!().\n
            Exiting interpolateLagrange!() without an interpolated solution."
    return nothing
  end

  if tf < 0.01
    @warn "tf needs to be greater than 0.01 to interpolate over solution.\n
          Exiting interpolateLagrange!() without an interpolated solution."
     return nothing
  end

  # sample points
  n.r.ocp.tpolyPts = [range(tf/n.ocp.Ni*(int-1),tf/n.ocp.Ni*int;length=numPts) .+ n.r.ocp.tst[1] for int in 1:n.ocp.Ni]
  n.r.ocp.XpolyPts = [[zeros(numPts) for int in 1:n.ocp.Ni] for st in 1:n.ocp.state.num]
  n.r.ocp.UpolyPts = [[zeros(numPts) for int in 1:n.ocp.Ni] for ctr in 1:n.ocp.control.num]
  if n.s.ocp.evalCostates; n.r.ocp.CSpolyPts = [[zeros(numPts) for int in 1:n.ocp.Ni] for st in 1:n.ocp.state.num] end

  # time data points
  t_st_int = [n.r.ocp.tst[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1] for int in 1:n.ocp.Ni]
  t_ctr_int = [n.r.ocp.tctr[n.ocp.Nck_cum[int]+1:n.ocp.Nck_cum[int+1]+1] for int in 1:n.ocp.Ni-1]
  t_ctr_int = push!(t_ctr_int, n.r.ocp.tst[n.ocp.Nck_cum[n.ocp.Ni]+1:n.ocp.Nck_cum[n.ocp.Ni+1]]) # -1 -> removing control in last mesh interval

  for int in 1:n.ocp.Ni
    # controls and states for this interval
    x_int, u_int = intervals(n, int, copy(n.r.ocp.X), copy(n.r.ocp.U))

     # sample polynomial in interval at n.r.ocp.tpolyPts
     for st in 1:n.ocp.state.num
      n.r.ocp.XpolyPts[st][int] = interpolate_lagrange(n.r.ocp.tpolyPts[int], t_st_int[int], x_int[:,st])'
     end

     for ctr in 1:n.ocp.control.num
      n.r.ocp.UpolyPts[ctr][int] = interpolate_lagrange(n.r.ocp.tpolyPts[int], t_ctr_int[int], u_int[:,ctr])'
     end

      # sample polynomial in interval at n.r.ocp.tpolyPts NOTE costate is missing the last point, that is the t_st_int[int][1:end-1]
      if n.s.ocp.evalCostates && n.s.ocp.evalConstraints
         for st in 1:n.ocp.state.num
          n.r.ocp.CSpolyPts[st][int] = interpolate_lagrange(n.r.ocp.tpolyPts[int], t_st_int[int][1:end-1], n.r.ocp.CS[st][int])'
         end
      end
  end

  # extract result into vectors
  temp = [n.r.ocp.tpolyPts[int][1:end] for int in 1:n.ocp.Ni] # time
  n.r.ocp.tpts = [idx for tempM in temp for idx=tempM]
  totalPts = length(n.r.ocp.tpts)

  n.r.ocp.Xpts = Matrix{Float64}(undef, totalPts, n.ocp.state.num)
  for st in 1:n.ocp.state.num # states
    temp = [n.r.ocp.XpolyPts[st][int][1:end,:] for int in 1:n.ocp.Ni]
    n.r.ocp.Xpts[:,st] = [idx for tempM in temp for idx=tempM]
  end

  n.r.ocp.Upts = Matrix{Float64}(undef, totalPts, n.ocp.control.num)
  for ctr in 1:n.ocp.control.num # controls
    temp = [n.r.ocp.UpolyPts[ctr][int][1:end,:] for int in 1:n.ocp.Ni]
    n.r.ocp.Upts[:,ctr] = [idx for tempM in temp for idx=tempM]
  end

  if n.s.ocp.evalCostates && n.s.ocp.evalConstraints
    n.r.ocp.CSpts = Matrix{Float64}(undef, totalPts, n.ocp.state.num)
    for st in 1:n.ocp.state.num # states
      temp = [n.r.ocp.CSpolyPts[st][int][1:end,:] for int in 1:n.ocp.Ni]
      n.r.ocp.CSpts[:,st] = [idx for tempM in temp for idx=tempM]
    end
  end

  return nothing
end

"""
"""
function interpolateLinear!(n::NLOpt; numPts::Int=250, tfOptimal::Any=false)
  # TODO throw an error if tfOptimal does not make sense given current solution
  if isa(tfOptimal,Bool)
    if n.s.ocp.finalTimeDV
      tf = getvalue(n.ocp.tf)
    else n.s.ocp.finalTimeDV
      tf = n.ocp.tf
    end
  else  # if there is a known optimal final time, then it is useful to evaluate the Lagrange polynomial at particular points to determine the error in the solution
    tf = tfOptimal
  end

  if isnan(tf)
      @warn "tf is a NaN cannot use it in interpolateLinear!().\n
            Exiting interpolateLinear!() without an interpolated solution."
      return nothing
  end

  if tf < 0.01
    @warn "tf needs to be greater than 0.01 to interpolate over solution.\n
          Exiting interpolateLinear!() without an interpolated solution."
    return nothing
  end

  # sample points
  t = range(0; length=numPts, stop=tf) .+ n.r.ocp.tst[1] # NOTE not tested
  n.r.ocp.tpts = convert(Array{Float64,1},t)
  n.r.ocp.Xpts = Matrix{Float64}(undef, numPts, n.ocp.state.num)
  n.r.ocp.Upts = Matrix{Float64}(undef, numPts, n.ocp.control.num)
  knots = (n.r.ocp.tst,)
  for st in 1:n.ocp.state.num
    sp_st = interpolate(knots,n.r.ocp.X[:,st],Gridded(Linear()))
    n.r.ocp.Xpts[:,st] = sp_st(n.r.ocp.tpts)
  end

  for ctr in 1:n.ocp.control.num
    if isequal(n.s.ocp.integrationMethod,:ps)
      sp_ctr = interpolate(knots,[n.r.ocp.U[:,ctr];0],Gridded(Linear()))
    else
      sp_ctr = interpolate(knots,n.r.ocp.U[:,ctr],Gridded(Linear()))
    end
    n.r.ocp.Upts[:,ctr] = sp_ctr(n.r.ocp.tpts)
  end
  return nothing
end

"""
plant2dfs!(n,sol)
"""
function plant2dfs!(n::NLOpt,sol,U)
  tSample = range(sol.t[1],sol.t[end],length=n.mpc.ip.state.pts)
  dfs = DataFrame()
  if isempty(n.r.ip.plant)
    dfs[:t] = tSample
    for st in 1:n.mpc.ip.state.num
      dfs[n.mpc.ip.state.name[st]] = [sol(t)[st] for t in tSample]
    end
    for ctr in 1:n.mpc.ip.control.num
      dfs[n.mpc.ip.control.name[ctr]] = [U[ctr][t] for t in tSample]
    end
  else
    dfs[!,:t] = [n.r.ip.plant[!,:t]; tSample]
    for st in 1:n.mpc.ip.state.num
      dfs[!,n.mpc.ip.state.name[st]] = [n.r.ip.plant[!,n.mpc.ip.state.name[st]]; [sol(t)[st] for t in tSample]]
    end
    for ctr in 1:n.mpc.ip.control.num
      dfs[!,n.mpc.ip.control.name[ctr]] = [n.r.ip.plant[!,n.mpc.ip.control.name[ctr]]; [U[ctr][t] for t in tSample] ]
    end
  end
  n.r.ip.plant = dfs

  # TODO only run this if saving, it is redundant data, or maybe put time stamps or something in the above data
  dfs = DataFrame()
  dfs[!,:t] = tSample
  for st in 1:n.mpc.ip.state.num
    dfs[!,n.mpc.ip.state.name[st]] = [sol(t)[st] for t in tSample]
  end
  for ctr in 1:n.mpc.ip.control.num
    dfs[!,n.mpc.ip.control.name[ctr]] = [U[ctr][t] for t in tSample]
  end
  push!(n.r.ip.dfsplant,dfs)
  return nothing
end

"""
dvs2dfs(n)
# funtionality to save state, costate, and control data from optimization
"""
function dvs2dfs(n::NLOpt)

    dfs = DataFrame()
    dfs[!, :t] = n.r.ocp.tst
    for st in 1:n.ocp.state.num
        dfs[!, n.ocp.state.name[st]] = n.r.ocp.X[:,st]
    end
    for ctr in 1:n.ocp.control.num
        if n.s.ocp.integrationMethod==:tm
            dfs[!, n.ocp.control.name[ctr]] = n.r.ocp.U[:,ctr]
        else
            dfs[!, n.ocp.control.name[ctr]] = [n.r.ocp.U[:,ctr];NaN]
        end
    end

    if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps && n.s.ocp.evalConstraints
        CS_vector = Matrix{Float64}(undef, n.ocp.state.pts, n.ocp.state.num)
        for st in 1:n.ocp.state.num # states
            temp = [n.r.ocp.CS[st][int][1:end,:] for int in 1:n.ocp.Ni]
            CS_vector[1:end-1,st] = [idx for tempM in temp for idx=tempM]
            CS_vector[end,st] = NaN
        end
        for st in 1:n.ocp.state.num # states
            dfs[!, Symbol(n.ocp.state.name[st],:cs)] = CS_vector[:,st]
        end
    end

    return dfs
end

"""
opt2dfs!(n)
# funtionality to save optimization data
"""
function opt2dfs!(n::NLOpt;kwargs...)

    kw = Dict(kwargs)

    if !haskey(kw,:statusUpdate)
        statusUpdate = false
    else
        statusUpdate = get(kw,:statusUpdate,0)
    end

    # make sure that the feildnames are initialized
    if isempty(n.r.ocp.dfsOpt)
        n.r.ocp.dfsOpt = DataFrame(tSolve = [], objVal = [], status = [], iterNum = [], evalNum = [])
    end

    if !statusUpdate
        push!(n.r.ocp.dfsOpt[!, :tSolve], n.r.ocp.tSolve)
        push!(n.r.ocp.dfsOpt[!, :objVal], n.r.ocp.objVal)
        push!(n.r.ocp.dfsOpt[!, :status], n.r.ocp.status)
    else  # TODO consider removing and cleaning this up
        push!(n.r.ocp.dfsOpt[!, :tSolve], NaN)
        push!(n.r.ocp.dfsOpt[!, :objVal], NaN)
        if statusUpdate && (typeof(n.r.ocp.status)==Symbol)
            push!(n.r.ocp.dfsOpt[!, :status], n.r.ocp.status)
        else
            push!(n.r.ocp.dfsOpt[!, :status], NaN)
        end
    end
    if n.s.mpc.on
        push!(n.r.ocp.dfsOpt[!, :iterNum], n.r.ocp.iterNum) # can set iter_nums for a higher/lower level algoritm
    end
    push!(n.r.ocp.dfsOpt[!, :evalNum], n.r.ocp.evalNum-1)

    return nothing
end

"""
con2dfs(n)
# funtionality to save constraint data
"""
function con2dfs(n::NLOpt)
  dfs_con=DataFrame()
  dfs_con[!, :conVal]=n.r.ocp.constraint.value
  return dfs_con
end

"""
postProcess!(n)
"""
function postProcess!(n::NLOpt; kwargs...)

    kw = Dict(kwargs)

    # check to see if the user is initializing while compensating for control delay
    if !haskey(kw,:Init)
        Init = false
    else
        Init = get(kw,:Init,0)
    end

    if n.s.ocp.save
        opt2dfs!(n)
    end

    # even if n.r.ocp.status==:Infeasible try to get solution. For the case that user may want to look at results to see where constraints where violated
    # in this case set =>  n.s.ocp.evalConstraints = true
    # http://jump.readthedocs.io/en/latest/refmodel.html#solve-status
    if !Init #&& (n.s.ocp.evalConstraints || ((n.r.ocp.status==:Optimal) || (n.r.ocp.status==:UserLimit)))
        if n.s.ocp.integrationMethod == :ps
            if n.s.ocp.finalTimeDV
                t = [scale_tau(n.ocp.ts[int],0.0,getvalue(n.ocp.tf)) for int in 1:n.ocp.Ni]     # scale time from [-1,1] to [t0,tf]
            else
                t = [scale_tau(n.ocp.ts[int],0.0,n.ocp.tf) for int in 1:n.ocp.Ni]
            end
            n.r.ocp.tctr = [idx for tempM in t for idx = tempM[1:end-1]] .+ getvalue(n.ocp.t0)
            n.r.ocp.tst = [n.r.ocp.tctr; t[end][end] .+ getvalue(n.ocp.t0)]
            # TODO check the above line... is t0 getting added on twice?
        elseif n.s.ocp.integrationMethod == :tm
            if n.s.ocp.finalTimeDV
                n.r.ocp.tctr = append!([0.0],cumsum(getvalue(n.ocp.dt))) .+ getvalue(n.ocp.t0)
            else
                n.r.ocp.tctr = append!([0.0],cumsum(n.ocp.dt)) .+ getvalue(n.ocp.t0)
            end
            n.r.ocp.tst = n.r.ocp.tctr
        end

        stateDataExists = false
        if n.r.ocp.status == :Optimal || (!n.s.mpc.onlyOptimal && n.s.mpc.on && isequal(n.mpc.v.evalNum,1))
            stateDataExists = true
            if (!isequal(n.r.ocp.status,:Optimal) && n.s.mpc.on && isequal(n.mpc.v.evalNum,1))
                @warn "There is no previous :Optimal solution to use since isequal(n.mpc.v.evalNum,1). \n
                    Attemting to extract: ",n.r.ocp.status," solution. \n
                    Setting: n.f.mpc.simFailed = [true, n.r.ocp.status] "
                n.f.mpc.simFailed = [true, n.r.ocp.status]
            end
            lengthStates = n.ocp.state.num * n.ocp.state.pts
            lengthControl = n.ocp.control.num * n.ocp.control.pts
            n.r.ocp.X = transpose(reshape(n.ocp.mdl.internalModel.inner.x[1:lengthStates], n.ocp.state.num, n.ocp.state.pts))
            n.r.ocp.U = transpose(reshape(n.ocp.mdl.internalModel.inner.x[lengthStates+1:lengthStates+lengthControl], n.ocp.control.num, n.ocp.control.pts))

        elseif n.s.mpc.on && n.s.mpc.lastOptimal && !n.s.mpc.onlyOptimal
            if !n.s.ocp.save
                error("This functionality currently needs to have n.s.ocp.save==true")
            end
            optIdx = findall(n.r.ocp.dfsOpt[:status].==:Optimal)[end]  # use the last :Optimal solution
            @show n.r.ocp.dfsOpt
            @show optIdx
            @show n.r.ocp.dfs[optIdx]
            @show n.mpc.v.t
            if n.r.ocp.dfs[optIdx][:t][1] > n.mpc.v.t
                timeIdx = 1
            elseif n.r.ocp.dfs[optIdx][:t][end] < n.mpc.v.t
                @warn "the current time is past the final time in the last :Optimal soultion.\n
                    Setting: n.f.mpc.simFailed = [true, n.r.ocp.status]"
                n.f.mpc.simFailed = [true, n.r.ocp.status]
                return nothing
            else
                timeIdx = findall(n.r.ocp.dfs[optIdx][:t] .- n.mpc.v.t .<= 0)[end]     # find the nearest index in time
            end
            # TODO: make an error message or fix      ERROR: LoadError: BoundsError: attempt to access 0-element Array{Int,1} at index [0]
            n.r.ocp.tst = n.r.ocp.dfs[optIdx][:t][timeIdx:end]
            n.r.ocp.X = zeros(Float64,length(n.r.ocp.dfs[optIdx][n.ocp.state.name[1]][timeIdx:end]),n.ocp.state.num)
            if n.s.ocp.integrationMethod==:tm  # TODO try to
                n.r.ocp.tctr = n.r.ocp.dfs[optIdx][:t][timeIdx:end]
                n.r.ocp.U = zeros(Float64,length(n.r.ocp.dfs[optIdx][n.ocp.control.name[1]][timeIdx:end]),n.ocp.control.num)
            else
                n.r.ocp.tctr = n.r.ocp.dfs[optIdx][:t][timeIdx:end-1]
                n.r.ocp.U = zeros(Float64,length(n.r.ocp.dfs[optIdx][n.ocp.control.name[1]][timeIdx:end-1]),n.ocp.control.num)
            end
            for st in 1:n.ocp.state.num
                n.r.ocp.X[:,st] = n.r.ocp.dfs[optIdx][n.ocp.state.name[st]][timeIdx:end]
            end
            for ctr in 1:n.ocp.control.num
                if n.s.ocp.integrationMethod==:tm
                n.r.ocp.U[:,ctr] = n.r.ocp.dfs[optIdx][n.ocp.control.name[ctr]][timeIdx:end]
                else
                n.r.ocp.U[:,ctr] = n.r.ocp.dfs[optIdx][n.ocp.control.name[ctr]][timeIdx:end-1]
                end
            end
        else
            @warn "The solution is not Optimal \n
                Setting: n.f.mpc.simFailed = [true, n.r.ocp.status] "
            n.f.mpc.simFailed = [true, n.r.ocp.status]
        end
        if n.s.ocp.evalConstraints && n.r.ocp.status!=:Error  # note may want to remove the && arg
            evalConstraints!(n)
            # TODO: make a note that costates can only be evaluated if .....
            if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps
                L1 = 0       # find index where dynamics constraints start
                for i in 1:length(n.r.ocp.constraint.name)
                    if n.r.ocp.constraint.name[i] == :dyn_con
                        L1 = n.r.ocp.constraint.nums[i][end][1]
                    end
                end
                mpb = JuMP.internalmodel(n.ocp.mdl)
                c = MathProgBase.getconstrduals(mpb)
                # NOTE for now since costates are not defined for :tm methods, n.r.ocp.CS is in a different format than n.r.ocp.X
                # in the future if costates are defined for :tm methods this can be changed
                n.r.ocp.CS = [[zeros(Float64,n.ocp.Nck[int]) for int in 1:n.ocp.Ni] for st in 1:n.ocp.state.num]
                for int in 1:n.ocp.Ni
                    b = 0
                    for st in 1:n.ocp.state.num
                        a = L1 + n.ocp.Nck[int]*(st-1)  # n.ocp.Nck[int]*(st-1) adds on indices for each additional state within the interval
                        b = a + n.ocp.Nck[int] - 1      # length of the state within this interval
                        n.r.ocp.CS[st][int] = -c[a:b]./n.ocp.ws[int]
                    end
                    L1 = b + 1 # adds indicies due to a change in the interval
                end
            end
        end

        if n.s.ocp.save && stateDataExists
            push!(n.r.ocp.dfs,dvs2dfs(n))
          # TODO figure out why the following line is broken
          #  push!(n.r.ocp.dfsCon,con2dfs(n))
            if n.s.ocp.interpolationOn
                if (n.r.ocp.status != :Error) # (n.r.ocp.status != :Infeasible) &&
                    if n.s.ocp.integrationMethod==:ps && (n.r.ocp.status != :Infeasible) && !n.s.ocp.linearInterpolation
                        interpolateLagrange!(n; numPts = n.s.ocp.numInterpPts, tfOptimal = n.s.ocp.tfOptimal)
                        push!(n.r.ocp.AlltpolyPts,n.r.ocp.tpolyPts)
                        push!(n.r.ocp.AllXpolyPts,n.r.ocp.XpolyPts)
                        push!(n.r.ocp.AllUpolyPts,n.r.ocp.UpolyPts)
                        if n.s.ocp.evalCostates && n.s.ocp.evalConstraints
                            push!(n.r.ocp.AllCSpolyPts,n.r.ocp.CSpolyPts)
                        end
                    else
                        interpolateLinear!(n; numPts = n.s.ocp.numInterpPts, tfOptimal = n.s.ocp.tfOptimal)
                        if n.s.ocp.integrationMethod==:ps && !n.s.ocp.linearInterpolation
                            push!(n.r.ocp.AlltpolyPts,nothing)
                            push!(n.r.ocp.AllXpolyPts,nothing)
                            push!(n.r.ocp.AllUpolyPts,nothing)
                            if n.s.ocp.evalCostates && n.s.ocp.evalConstraints
                                push!(n.r.ocp.AllCSpolyPts,nothing)
                            end
                        end
                    end
                end
            end  # TODO: save [] of interpolated pts for :tm methods
        end
    end

    return nothing
end
# TODO: add some "finalPostProcess" to save time

"""
optimize!(n)
# solves JuMP model and saves optimization data
"""
function optimize!(n::NLOpt; Iter::Int=0)

    t1 = time()
    status = JuMP.solve(n.ocp.mdl)
    t2 = time()

    if !n.s.ocp.cacheOnly
        n.r.ocp.status = status
        n.r.ocp.tSolve = t2 - t1
        n.r.ocp.objVal = getobjectivevalue(n.ocp.mdl)
        n.r.ocp.iterNum = Iter    # possible iteration number for a higher level algorithm
        n.r.ocp.evalNum = n.r.ocp.evalNum + 1
        postProcess!(n)      # temporarily save data
    end

    return nothing
end

"""
"""
function evalConstraints!(n::NLOpt)

    n.r.ocp.constraint.value = []   # reset values
    n.r.ocp.constraint.nums = []
    s=1

    for i = 1:length(n.r.ocp.constraint.handle)
        if n.r.ocp.constraint.name[i] == :dyn_con  # state constraits
            dfs=Vector{DataFrame}(n.ocp.state.num)
            con=DataFrame(step=1)
            l=0
            for st in 1:n.ocp.state.num
                if n.s.ocp.integrationMethod==:ps
                    temp= [ getdual(n.r.ocp.constraint.handle[i][int][:,st]) for int in 1:n.ocp.Ni ]
                    vals = [ idx for tempM in temp for idx in tempM ]
                    dfs[st] = DataFrame(step=1:sum(n.ocp.Nck) ; Dict(n.ocp.state.name[st] => vals)...)
                    l += length(vals)
                else
                    dfs[st] = DataFrame(step=1:n.ocp.N; Dict(n.ocp.state.name[st] => getdual(n.r.ocp.constraint.handle[i][:,st]))...)
                    l += length(n.r.ocp.constraint.handle[i][:,st])
                end

                con = ( st == 1 ? dfs[st] : join(con, dfs[st], on=:step) )

            end
        else
            S=0
            try
                S=JuMP.size(n.r.ocp.constraint.handle[i])
            catch
                error("\n For now, the constraints cannot be in this form: \n
                con=@NLconstraint(mdl,n.r.ocp.u[1,1]==param); \n
                Write it in array form: \n
                con=@NLconstraint(mdl,[i=1],n.r.ocp.u[i,1]==param); \n")
            end
            if length(S) == 1
                con = DataFrame(step=1:length(n.r.ocp.constraint.handle[i]); Dict(n.r.ocp.constraint.name[i] => getdual(n.r.ocp.constraint.handle[i][:]))...)
                l=S[1]
            elseif length(S) == 2
                dfs=Vector{DataFrame}(S[1])
                con=DataFrame(step=1)
                for idx in 1:S[1]
                    try
                        dfs[idx] = DataFrame(step=1:S[2]; Dict(n.r.ocp.constraint.name[i] => getdual(n.r.ocp.constraint.handle[i][idx,:]))...)
                    catch
                        dfs[idx] = DataFrame(step=1:S[2]; Dict(n.r.ocp.constraint.name[i] => NaN)...) # fix for case where all of the states are not being constrainted, but some are within some XF_tol
                    end
                    if idx==1;
                        con = dfs[idx]
                    else
                        con = join(con,dfs[idx],on=:step,makeunique=true)
                    end
                end

                l = S[1] * S[2]
            end
        end
        f = s + l - 1
        num = (i, n.r.ocp.constraint.name[i], "length = $l", string("indices in g(x) = "), (s, f))
        push!(n.r.ocp.constraint.nums,num)
        push!(n.r.ocp.constraint.value,con)
        s = f + 1
    end
    return nothing
end

"""
description = string(
" *  \n ")

resultsDir!(r;resultsName=resultsName,description=description)
# removes results folder and creates a new one
# TODO consider putting in a warning or some sort of an interaction with user
"""
function resultsDir!(n; resultsName::String = "",description::DataFrame = DataFrame())

    results_dir=string(n.r.mainDir,"/results/",resultsName)  # define directories
    n.r.resultsDir = results_dir

    if isdir(n.r.resultsDir)
        rm(n.r.resultsDir; recursive=true)
        print("\n The old results have all been deleted! \n \n")
    end
    mkdir(n.r.resultsDir) # create directory

    cd(n.r.resultsDir)
    CSV.write("description.csv", description; quotechar = ' ')
    cd(n.r.mainDir)
    return nothing
end

########################################################################################
# MPC functions
########################################################################################
"""
defineMPC!(n)
"""
function defineMPC!(n::NLOpt;
                   mode::Symbol=:OCP,
                   predictX0::Bool=true,
                   fixedTp::Bool=true,
                   tp::Any=Any,
                   tex::Float64=0.5,
                   IPKnown::Bool=true,
                   saveMode::Symbol=:all,
                   maxSim::Int=100,
                   goal=n.ocp.XF,
                   goalTol= 0.1 * abs.(n.ocp.X0 - n.ocp.XF), # TODO: Set default goal tolerance elsewhere?
                   lastOptimal::Bool=true,
                   printLevel::Int=2,
                   onlyOptimal::Bool=false)

    n.s.mpc.on = true

    #n.mpc::MPC{T} = MPC()
    n.s.mpc.mode = mode
    n.s.mpc.predictX0 = predictX0
    n.s.mpc.fixedTp = fixedTp
    n.mpc.v.tp = tp
    n.mpc.v.tex = tex
    n.s.mpc.IPKnown = IPKnown
    n.s.mpc.saveMode = saveMode
    n.s.mpc.maxSim = maxSim
    n.mpc.v.goal = goal
    n.mpc.v.goalTol = goalTol
    n.s.mpc.lastOptimal = lastOptimal
    n.s.mpc.printLevel = printLevel
    n.s.mpc.onlyOptimal = onlyOptimal

    n.f.mpc.simFailed[1] = false # TODO: figure out why this is getting defined as 0.0, not false during initialization
    n.f.mpc.defined = true

    return nothing
end

"""
# TODO consider letting user pass options
"""
function initOpt!(n::NLOpt; save::Bool=true, evalConstraints::Bool=false)

    if n.s.mpc.on
        error("call initOpt!() before defineMPC!(). initOpt!() will destroy n")
    end

    n.s.ocp.save = false
    n.s.mpc.on = false
    n.s.ocp.evalConstraints = false
    n.s.ocp.cacheOnly = true

    if n.s.ocp.save
        @warn "saving initial optimization results where functions where cached!"
    end

    for k in 1:n.mpc.v.initOptNum # initial optimization (s)
        status = optimize!(n)
        if status == :Optimal
            break;
        end
    end

    # defineSolver!(n,solverConfig(c)) # modifying solver settings NOTE currently not in use

    n.s.ocp.save = save  # set to false if running in parallel to save time
    n.s.ocp.cacheOnly = false
    n.s.ocp.evalConstraints = evalConstraints # set to true to investigate infeasibilities

    return nothing
end

"""
# add a mode that solves as quickly as possible
# consider using the IP always.
defineModel!(n)
"""
function defineIP!(n::NLOpt,model;stateNames=[],controlNames=[],X0a=[])

   if n.s.mpc.mode == :OCP # this function is called automatically for this mode
    if !isempty(stateNames)
        error("stateNames are set automatically for :mode == :OCP and cannot be provided.")
    end
    if !isempty(controlNames)
        error("controlNames are set automatically for :mode == :OCP and cannot be provided.")
    end
    if !isempty(X0a)
        error("X0a is set automatically for :mode == :OCP and cannot be provided.")
    end
    n.r.ip.X0a = copy(n.ocp.X0)  # NEED to append time
    n.mpc.ip.state.model = model
    n.mpc.ip.state.name = n.ocp.state.name
    n.mpc.ip.state.description = n.ocp.state.description
    n.mpc.ip.state.num = n.ocp.state.num
    n.mpc.ip.state.pts = n.ocp.state.pts

    n.mpc.ip.control.name = n.ocp.control.name
    n.mpc.ip.control.description = n.ocp.control.description
    n.mpc.ip.control.num = n.ocp.control.num
    n.mpc.ip.control.pts = n.ocp.control.pts

    # add X0 t0 plant dfs
    n.r.ip.plant[!, :t] = [n.mpc.v.t0]
    for st in 1:n.mpc.ip.state.num
      n.r.ip.plant[!, n.mpc.ip.state.name[st]] = [copy(n.ocp.X0)[st]]
    end
    for ctr in 1:n.mpc.ip.control.num
      n.r.ip.plant[!, n.mpc.ip.control.name[ctr]] = [0]
    end

   elseif isequal(n.s.mpc.mode,:IP)
    if isempty(stateNames)
     error("unless :mode == :OCP the stateNames must be provided.")
    end
    if isempty(controlNames)
     error("unless :mode == :OCP the controlNames must be provided.")
    end
    if isempty(X0a)
     error("unless :mode == :OCP X0a must be provided.")
    end
    if isempty(model)
     error("A model needs to be passed for the IP mode.")
    else
    if isequal(length(X0a),length(stateNames))
      error(string("\n Length of X0a must match length(stateNames) \n"))
    end

     n.mpc.ip.state::State = State() # reset
     n.mpc.ip.state.num = length(stateNames)
     for i in 1:n.mpc.ip.state.num
       if stateNames[i]==:xxx
         error("xxx is OFF limits for a state name; please choose something else. \n")
       end
       push!(n.mpc.ip.state.name,stateNames[i])
     end

     n.mpc.ip.control::Control = Control() # reset
     n.mpc.ip.control.num = length(controlNames)
     for i in 1:n.mpc.ip.control.num
       if controlNames[i]==:xxx
         error("xxx is OFF limits for a control name; please choose something else. \n")
       end
       push!(n.mpc.ip.control.name,controlNames[i])
     end
     n.mpc.r.ip.X0a = X0a
     n.mpc.ip.state.model = model # TODO validate typeof model
    end
   elseif isequal(n.s.mpc.mode,:EP)
    error("not setup for :EP")
   else
    error("n.mpc.s.mode = ",n.s.mpc.mode," not defined." )
   end

   # consider calling mapNames
   return nothing
end

"""
mapNames!(n)
"""
function mapNames!(n::NLOpt)
  if isequal(n.s.mpc.mode,:IP)
    s1 = n.ocp.state.name
    c1 = n.ocp.control.name
    s2 = n.mpc.ip.state.name
    c2 = n.mpc.ip.control.name
  elseif isequal(n.s.mpc.mode,:EP)
    error(":EP function not ready")
  else
    error("mode must be either :IP or :EP")
  end

  m = []
  # go through all states in OCP
  idxOCP = 1
  for var in s1
    # go through all states in IP
    idxIP = findall(var.==s2)
    if !isempty(idxIP)
      push!(m, [var; :stOCP; idxOCP; :stIP; idxIP[1]])
    end

    # go through all controls in IP
    idxIP = findall(var.==c2)
    if !isempty(idxIP)
      push!(m, [var; :stOCP; idxOCP; :ctrIP; idxIP[1]])
    end
    idxOCP = idxOCP + 1
  end

  # go through all controls in OCP
  idxOCP = 1
  for var in c1
    # go through all states in IP
    idxIP = findall(var.==s2)
    if !isempty(idxIP)
      push!(m, [var; :ctrOCP; idxOCP; :stIP; idxIP[1]])
    end

    # go through all controls in IP
    idxIP = findall(var.==c2)
    if !isempty(idxIP)
      push!(m, [var; :ctrOCP; idxOCP; :ctrIP; idxIP[1]])
    end
    idxOCP = idxOCP + 1
  end

  if isequal(n.s.mpc.mode,:IP)
    n.mpc.mIP = m
  elseif isequal(n.s.mpc.mode,:EP)
    error(":EP function not ready")
  else
    error("mode must be either :IP or :EP")
  end

  return nothing
end

"""
# TODO fix this so that prediction simulation time is ahead. May not effect results.
# NOTE this may be ok... we are getting X0p for initialization of the OCP
# as long as the OCP pushes the time ahead. which it does then everything is fine!!
# consider making user pass X0, t0, tf
"""
function simIPlant!(n::NLOpt)
  if isequal(n.mpc.ip.state.pts,0)
   error("isqual(n.mpc.ip.state.pts,0), cannot simulate with zero points.")
  end
  X0 = currentIPState(n)[1]
  t0 = round(n.mpc.v.t, digits=3) # if rounding is too rough, then the tex will be essentially 0!
  tf = round(n.mpc.v.t + n.mpc.v.tex, digits=3)

  if isequal(n.s.mpc.mode,:OCP)
   if isequal(n.mpc.v.evalNum,1)
    U = 0*Matrix{Float64}(undef, n.ocp.control.pts,n.ocp.control.num)
    t = Vector(range(t0,tf,length=n.ocp.control.pts))
   elseif n.s.ocp.interpolationOn
    U = n.r.ocp.Upts
    t = n.r.ocp.tpts
   else
    U = n.r.ocp.U
    t = n.r.ocp.tctr
   end
  else
   error("TODO")
  end
  # chop of first control point for bkwEuler as it is typically 0
  if isequal(n.s.ocp.integrationScheme,:bkwEuler)
   U = U[2:end,:]
   t = t[2:end]
  end

  sol, U = n.mpc.ip.state.model(n,X0,t,U,t0,tf)
  return sol, U
end

"""
"""
function currentIPState(n::NLOpt)
  if isempty(n.r.ip.plant)
    error("there is no data in n.r.ip.plant")
  end

  # even though may have solution for plant ahead of time
  # can only get the state up to n.mpc.v.t
  idx = findall((n.mpc.v.t .- n.r.ip.plant[!,:t]) .>= 0)
  if isempty(idx)
    error("(n.mpc.v.t - n.r.ip.plant[!,:t]) .>= 0) is empty.")
  else
    X0 = [zeros(n.mpc.ip.state.num),n.mpc.v.t]
    for st in 1:n.mpc.ip.state.num
      X0[1][st] = n.r.ip.plant[!, n.mpc.ip.state.name[st]][idx[end]]
    end
  end
  return X0
end

"""
# TODO eventually the "plant" will be different from the "model"
predictX0!(n)
"""
function predictX0!(n::NLOpt)

  if n.s.mpc.fixedTp
   # NOTE consider passing back (n.mpc.v.t + n.mpc.v.tex) from simIPlant!()
   tp = round(n.mpc.v.t + n.mpc.v.tex, digits=1)  # TODO add 1 as an MPCparamss
  else
   error("TODO")
  end

  if isequal(n.s.mpc.mode,:OCP)
   sol, U = simIPlant!(n)
   X0p = [sol(sol.t[end])[:],tp]
   push!(n.r.ip.X0p,X0p)
  else
    error("TODO")
  end
 # else
   # with no control signals to follow, X0p is simply the current known location of the plant
 #  X0 = currentIPState(n)
 #  X0p = [X0[1], tp]  # modify X0 to predict the time
  # push!(n.r.ip.X0p,X0p)
 # end
  return nothing
end

"""
"""
function updateX0!(n::NLOpt,args...)
 # need to map n.r.ip.X0p to n.X0 (states may be different)
 # NOTE for the :OCP mode this is OK

 if !n.s.mpc.predictX0 #  use the current known plant state to update OCP
   push!(n.r.ip.X0p,currentIPState(n))
 else
   predictX0!(n)
 end

 if isequal(n.s.mpc.mode,:OCP)
   if !isequal(length(args),0)
    X0 = args[1]
    if length(X0)!=n.ocp.state.num
      error(string("\n Length of X0 must match number of states \n"));
    end
    n.ocp.X0 = X0
   else
    n.ocp.X0 = n.r.ip.X0p[end][1] # the n.ocp. structure is for running things
   end
   push!(n.r.ocp.X0, n.ocp.X0)    # NOTE this may be for saving data
   setvalue(n.ocp.t0, copy(n.r.ip.X0p[end][2]))
 else
  error("not set up for this mode")
 end

  if n.s.mpc.shiftX0 # TODO consider saving linear shifting occurances
    for st in 1:n.ocp.state.num
      if n.ocp.X0[st] < n.ocp.XL[st]
        n.ocp.X0[st] = n.ocp.XL[st]
      end
      if n.ocp.X0[st] > n.ocp.XU[st]
        n.ocp.X0[st] = n.ocp.XU[st]
      end
    end
  end
  # update states with n.ocp.X0
  for st in 1:n.ocp.state.num
    if n.s.ocp.x0slackVariables
     JuMP.setRHS(n.r.ocp.x0Con[st,1], n.ocp.X0[st])
     JuMP.setRHS(n.r.ocp.x0Con[st,2],-n.ocp.X0[st])
    else
      JuMP.setRHS(n.r.ocp.x0Con[st],n.ocp.X0[st])
    end
  end
  return nothing
end


"""
"""
function goalReached!(n::NLOpt,args...)
  if isequal(n.s.mpc.mode,:OCP)
    X = currentIPState(n)[1]
  else
    X = args[1]
    error("TODO")
  end
  A = (abs.(X - n.mpc.v.goal) .<= n.mpc.v.goalTol)
  B = isnan.(n.mpc.v.goal)
  C = [A[i]||B[i] for i in 1:length(A)]

  if all(C)
   if isequal(n.s.mpc.printLevel,2)
    println("Goal Attained! \n")
   end
    n.f.mpc.goalReached = true
  elseif n.s.mpc.expandGoal && (getvalue(n.ocp.tf) < n.mpc.v.tex)
    A =( abs.(X - n.mpc.v.goal) .<= n.s.mpc.enlargeGoalTolFactor*n.mpc.v.goalTol)
    C = [A[i]||B[i] for i in 1:length(A)]
    if all(C)
     if isequal(n.s.mpc.printLevel,2)
      println("Expanded Goal Attained! \n")
     end
     n.f.mpc.goalReached = true
    else
     println("Expanded Goal Not Attained! \n
              Stopping Simulation")
     n.f.mpc.simFailed = [true, :expandedGoal]
    end
  end

 return n.f.mpc.goalReached
end
# if the vehicle is very close to the goal sometimes the optimization returns with a small final time
# and it can even be negative (due to tolerances in NLP solver). If this is the case, the goal is slightly
# expanded from the previous check and one final check is performed otherwise the run is failed
#if getvalue(n.ocp.tf) < 0.01
#  if ((n.r.ip.dfplant[end][:x][end]-c["goal"]["x"])^2 + (n.r..ip.dfplant[end][:y][end]-c["goal"]["yVal"])^2)^0.5 < 2*c["goal"]["tol"]
#  println("Expanded Goal Attained! \n"); n.f.mpc.goal_reached=true;
#  break;/
#  else
#  warn("Expanded Goal Not Attained! -> stopping simulation! \n"); break;
#  end
#elseif getvalue(n.ocp.tf) < 0.5 # if the vehicle is near the goal => tf may be less then 0.5 s
#  tf = (n.r.evalNum-1)*n.mpc.v.tex + getvalue(n.ocp.tf)
#else
#  tf = (n.r.evalNum)*n.mpc.v.tex
#end


"""
"""
function simMPC!(n::NLOpt;updateFunction::Any=[],checkFunction::Any=[])
  for ii = 1:n.s.mpc.maxSim
    if isequal(n.s.mpc.printLevel,2)
     println("Running model for the: ",n.mpc.v.evalNum," time")
    end
    #############################
    # (A) and (B) in "parallel"
    #############################

    # (B) simulate plant
    sol, U = simIPlant!(n) # the plant simulation time will lead the actual time
    plant2dfs!(n,sol,U)

    # check to see if the simulation failed (i.e. plant crashed)
    if !isequal(typeof(checkFunction),Array{Any,1})
     val, sym = checkFunction(n)
     if val
       n.f.mpc.simFailed = [val, sym]
      break
     end
    end

    # check to see if the goal has been reached
    if goalReached!(n); break; end
    if n.f.mpc.simFailed[1]; break; end

    # (A) solve OCP  TODO the time should be ahead here as it runs
    updateX0!(n)  # before updateFunction()
    if !isequal(typeof(updateFunction),Array{Any,1})
      updateFunction(n)
    end

    optimize!(n)
    if n.f.mpc.simFailed[1]; break; end

    # advance time
    n.mpc.v.t = n.mpc.v.t + n.mpc.v.tex
    n.mpc.v.evalNum = n.mpc.v.evalNum + 1
  end
end
# practical concerns/questions
#################################
# 1) predict X0
# 2) shift X0 for NLP feasibility
# 3) ensuring that U passed to the plant is feasible
     # effected by interpolation, demonstrate by varying numPts
     # seems to be a major problem with LGR nodes, possibly due to Runge effect
# 4) fixedTp or variableTp
# 5) usePrevious optimal.
     # at some point will be unable to do this
# 6) infeasibilities, soft constraint on inital conditions

# TODO
# 1) plot the goal, the tolerances on X0p
# 2) calculate the error and plot


"""
linearStateTolerances!(n::NLOpt)
# the purpose of this function is to taper the tolerances on the constant state constraints
# the idea is that when doing MPC, the final states are well within the bounds so that the next optimization is not initalized at an infeasible point
# if you want a constant bond, set the slope to zero
# default is a positive slope on the lower bound and a negative slope on the upper bound
# this functionality in not needed for states like position, so you do not need to add a linearStateTol for all states

"""
function linearStateTolerances!(n::NLOpt;
                                mXL::Array{Any,1}=falses(n.ocp.state.num),
                                mXU::Array{Any,1}=falses(n.ocp.state.num))
  n.ocp.mXL=mXL;  n.ocp.mXU=mXU;
  for st in 1:n.ocp.state.num
    # lower state constraint
    if n.ocp.mXL[st]!=false
      if !isnan(n.ocp.XL[st])
        for j in 1:n.ocp.state.pts
            n.ocp.XL_var[st,j]=n.ocp.XL[st] + n.ocp.mXL[st]*(j/n.ocp.state.pts); # lower
        end
      end
    end

    # upper state constraint
    if n.ocp.XU[st]!=false
      if !isnan(n.ocp.XU[st])
        for j in 1:n.ocp.state.pts
        n.ocp.XU_var[st,j]=n.ocp.XU[st] + n.ocp.mXU[st]*(j/n.ocp.state.pts); # upper
        end
      end
    end
  end
  return nothing
end

"""
defineTolerances!(n::NLOpt)
"""
function defineTolerances!(n::NLOpt;
                          X0_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,),
                          XF_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,))
  # TODO error checking, if the user does not pass tolerances etc.
  n.ocp.X0_tol = X0_tol
  n.ocp.XF_tol = XF_tol
  return nothing
end

"""
create_tV!(n::NLOpt)
# define a time vector (n.ocp.tV) for use with time varying constraints when (finalTimeDV=>true)
"""
function create_tV!(n::NLOpt)

  if n.s.ocp.integrationMethod==:ps
    # create mesh points, interval size = tf_var/Ni
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.Ni+1], (idx-1)*n.ocp.tf/n.ocp.Ni)
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts = [Array{Any}(undef, n.ocp.Nck[int]+1,) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      ts[int][1:end-1] = @NLexpression(n.ocp.mdl,[j=1:n.ocp.Nck[int]], (tm[int+1]-tm[int])/2*n.ocp.tau[int][j] +  (tm[int+1]+tm[int])/2);
      ts[int][end] = @NLexpression(n.ocp.mdl, n.ocp.tf/n.ocp.Ni*int) # append +1 at end of each interval
    end
    tt1 = [idx for tempM in ts for idx = tempM[1:end-1]]
    tmp = [tt1;ts[end][end]]
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
  else
    # create vector with the design variable in it
    t = Array{Any}(undef, n.ocp.N+1,1)
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.N], n.ocp.tf/n.ocp.N*idx)
    tmp = [0;tm]
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.ocp.t0 + tmp[j])
  end
  return nothing
end

"""
obj=integrate!(n::NLOpt,:(u1))
"""
function integrate!(n::NLOpt,V::Expr)
  if n.s.ocp.integrationMethod==:ps
    integral_expr = [Array{Any}(undef, n.ocp.Nck[int]) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)
      L = size(x_int)[1]-1
      integral_expr[int][:] = NLExpr(n,V,x_int,u_int,L)
    end
    @NLexpression(n.ocp.mdl, temp[int=1:n.ocp.Ni], (n.ocp.tf-n.ocp.t0)/2*sum(n.ocp.ws[int][j]*integral_expr[int][j] for j = 1:n.ocp.Nck[int]) )
    expression = @NLexpression(n.ocp.mdl, sum(temp[int] for int = 1:n.ocp.Ni))
  elseif n.s.ocp.integrationMethod==:tm
    L = size(n.r.ocp.x)[1]
    temp = NLExpr(n,V,n.r.ocp.x,n.r.ocp.u,L);
    if n.s.ocp.integrationScheme==:bkwEuler
      # NOTE integration this way does not penalize the first control
      expression = @NLexpression(n.ocp.mdl, sum(temp[j+1]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
      #expression = @NLexpression(n.ocp.mdl, sum(temp[j]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    elseif n.s.ocp.integrationScheme==:trapezoidal
      expression = @NLexpression(n.ocp.mdl, sum(0.5*(temp[j]+temp[j+1])*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    else
      error("\n $(n.s.ocp.integrationScheme) not defined in integrationSchemes\n")
    end
  else
    error("\n $(n.s.ocp.integrationMethod) not defined in integrationMethods \n")
  end
  return expression
end

########################################################################################
# constraint data functions
########################################################################################
"""
# funtionality to save constraint data
"""
function initConstraint!(n::NLOpt)
    if n.r.ocp.constraint === nothing
        n.r.ocp.constraint = Constraint()
    end
    return nothing
end

"""
"""
function newConstraint!(n::NLOpt, handle, name::Symbol)
    initConstraint!(n)
    # @info "\nEntering newConstraint!\n"
    # @info "newConstraint!: size of n.r.ocp.constraint.handle:   $(size(n.r.ocp.constraint.handle))"
    # @info "newConstraint!: type of n.r.ocp.constraint.handle:   $(typeof(n.r.ocp.constraint.handle))"
    # @info "newConstraint!: size of n.r.ocp.constraint.name:     $(size(n.r.ocp.constraint.name))"
    # @info "newConstraint!: type of n.r.ocp.constraint.name:     $(typeof(n.r.ocp.constraint.name))"
    # @info "newConstraint!: size of handle:                      $(size(handle))"
    # @info "newConstraint!: type of handle:                      $(typeof(handle))"
    # @info "newConstraint!: value of handle:                     $handle"
    #@info "newConstraint!: type of name:                        $(typeof(name))"
    # @info "newConstraint!: value of name:                       $name"
    # TODO test if this is the test that is making sure constraints can be entered in results in postProcess
    #@info "newConstraint!: testing vcat:                        $(vcat(n.r.ocp.constraint.handle, handle))"

    push!(n.r.ocp.constraint.name,name)
    #n.r.ocp.constraint.handle = vcat(n.r.ocp.constraint.handle, handle)
    push!(n.r.ocp.constraint.handle,handle)

    return nothing
end

"""
maxDualInf = evalMaxDualInf(n::NLOpt)
# funtionality to evaluate maximum dual infeasibility of problem
"""
function evalMaxDualInf(n::NLOpt)
    num=length(n.r.ocp.constraint.handle); dual_con_temp=zeros(num);
    for i = 1:length(n.r.ocp.constraint.handle)
        if !isempty(n.r.ocp.constraint.handle[i])
            if n.r.ocp.constraint.name[i]==:dyn_con  # state constraits
                temp1=zeros(n.ocp.state.num);
                for st in 1:n.ocp.state.num
                    if n.s.ocp.integrationMethod==:ps
                        temp=[getdual(n.r.ocp.constraint.handle[i][int][:,st]) for int in 1:n.ocp.Ni];
                        vals=[idx for tempM in temp for idx=tempM];
                        temp1[st]=maximum(vals);
                    else
                        temp1[st] = maximum(getdual(n.r.ocp.constraint.handle[i][:,st]));
                    end
                end
                dual_con_temp[i]=maximum(temp1);
            else
                S=JuMP.size(n.r.ocp.constraint.handle[i])
                if length(S)==1
                dual_con_temp[i]=maximum(getdual(n.r.ocp.constraint.handle[i][:]));
                elseif length(S)==2
                temp1=zeros(S[1]);
                for idx in 1:S[1]
                    temp1[idx] = maximum(getdual(n.r.ocp.constraint.handle[i][idx,:]));
                end
                dual_con_temp[i]=maximum(temp1);
                end
            end
        end
    end
    return maximum(maximum(maximum(dual_con_temp)))
end

########################################################################################
# state data functions
########################################################################################

"""
# initialize states
"""
function initState(numStates)::State
  s = State()
  s.num = numStates
  s.name = [Symbol("x$i") for i in 1:numStates]
  s.description = [String("x$i") for i in 1:numStates]
  return s
end

"""
"""
function states!(n::NLOpt,names;descriptions=[])
  if !n.f.ocp.defined
    error("\n call define() before calling states!() \n")
  end
  if length(names)!=n.ocp.state.num
    error("\n Check size of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.state.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.state.num
    if names[i]==:xxx
      error("xxx is OFF limits for a state name; please choose something else. \n")
    end
  end

   n.ocp.state.name = names
   if !isempty(descriptions)
     n.ocp.state.description = descriptions
   end
  return nothing
end

########################################################################################
# control data functions
########################################################################################
"""
# initialize control
"""
function initControl(numControls)::Control
  c = Control()
  c.num = numControls
  c.name = [Symbol("u$i") for i in 1:numControls]
  c.description = [String("u$i") for i in 1:numControls]
  return c
end

"""
"""
function controls!(n::NLOpt,names;descriptions=[])
  if !n.f.ocp.defined
    error("\n call define() before calling controls!() \n")
  end
  if length(names)!=n.ocp.control.num
    error("\n Check sizes of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.control.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.control.num
    if names[i]==:uuu
      error("uuu is OFF limits for a control name; please choose something else. \n")
    end
   end

  n.ocp.control.name = names
  if !isempty(descriptions)
    n.ocp.control.description = descriptions
  end
  return nothing
end

########################################################################################
# data related functions
########################################################################################

"""
"""
function savePlantData!(n::NLOpt)
  if isempty(n.r.ip.plant)
    @warn "No plant data to save.\n
          Make sure that n.s.save = true "
    return nothing
  end
  cd(n.r.resultsDir)
    CSV.write("plant.csv", n.r.ip.plant; quotechar = ' ');
  cd(n.r.mainDir)
  return nothing
end

"""
"""
function saveOptData(n::NLOpt)
  if isempty(n.r.ocp.dfsOpt)
    @warn "No optimization data to save.\n
          Make sure that n.s.save = true "
    return nothing
  end
  cd(n.r.resultsDir)
    CSV.write("opt.csv", n.r.ocp.dfsOpt; quotechar = ' ');
  cd(n.r.mainDir)
  return nothing
end


"""
"""
function saveData(n::NLOpt)
  # all polynomial data
  dfs=DataFrame();

  dfs[!,:t] = n.r.ocp.tpts
  for st in 1:n.ocp.state.num # state
    dfs[n.ocp.state.name[st]]=n.r.ocp.Xpts[:,st];
  end

  for ctr in 1:n.ocp.control.num # control
    dfs[n.ocp.control.name[ctr]]=n.r.ocp.Upts[:,ctr];
  end

  if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps && n.s.ocp.evalConstraints
    for st in 1:n.ocp.state.num # state
      dfs[Symbol(n.ocp.state.name[st],:cs)]=n.r.ocp.CSpts[:,st];
    end
  end

  cd(n.r.resultsDir)
    CSV.write("st_ctr.csv",n.r.ocp.dfs[end]; quotechar = ' '); # assuming only want the last one is needed
    CSV.write("st_ctr_poly.csv",dfs; quotechar = ' ');
  cd(n.r.mainDir)
  return nothing
end
"""
"""
function saveBenchMarkData!(n::NLOpt)
  first=2
  dfs=DataFrame();
  temp = [n.r.ocp.dfs[jj][:t][1:end-1,:] for jj in first:length(n.r.ocp.dfs)]; # time
  U=[idx for tempM in temp for idx=tempM]; dfs[:t]=U;

  for st in 1:n.ocp.state.num # state
    temp = [n.r.ocp.dfs[jj][n.ocp.state.name[st]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.ocp.state.name[st]]=U;
  end

  for ctr in 1:n.ocp.control.num # control
    temp = [n.r.ocp.dfs[jj][n.ocp.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.ocp.control.name[ctr]]=U;
  end

  # save optimization times
  temp = [n.r.ocp.dfsOpt[jj][n.ocp.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];

  cd(n.r.resultsDir)
    CSV.write("bench_data.csv",dfs; quotechar = ' ');
  cd(n.r.mainDir)
  return nothing
end


"""
# maximum(x->maximum(x[:A]), dfs) -> consider
# maximum(x->maximum(filter(y->y !== nothing, x[:A])), dfs)
# minimum(x->maximum(filter(y->y !== nothing, x[:t])), dfs)
"""
function minDF(dfs,varb)
  k=length(dfs); tmp=1e10*ones(k,1);
  for i in 1:k
    if dfs[i] !== nothing
      tmp[i]=minimum(dfs[i][varb])
    end
  end
  minimum(tmp)  # find the minimum
end

"""
"""
function maxDF(dfs,varb)
  k=length(dfs); tmp=1e-10*ones(k,1);
  for i in 1:k
    if dfs[i] !== nothing
      tmp[i]=maximum(dfs[i][varb])
    end
  end
  maximum(tmp)  # find the minimum
end

"""
"""
function try_import(name::Symbol)
    try
        @eval using $name
        return true
    catch e
        error("Could not import \"$name\" into Julia environment.\nOriginal Error:\n$e")
    end
end
