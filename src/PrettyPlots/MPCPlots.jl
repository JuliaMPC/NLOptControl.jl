"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 04/26/2018 Last Modified: 04/26/2018 \n
--------------------------------------------------------------------------------------\n
"""
function mpcPlot(n,idx,kwargs...)
  if !isdir(n.r.resultsDir); resultsDir!(n); end
  stp = [statePlot(n,idx,st;kwargs...) for st in 1:n.ocp.state.num]
  ctp = [controlPlot(n,idx,ctr;kwargs...) for ctr in 1:n.ocp.control.num]

  tp = tPlot(n,idx)

  if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps && n.s.ocp.evalConstraints
    csp = [costatePlot(n,idx,st;kwargs...) for st in 1:n.ocp.state.num]
    all = [stp;ctp;tp;csp]
  else
    all = [stp;ctp;tp]
  end

  h = plot(all...,size=_pretty_defaults[:size])
  if !_pretty_defaults[:simulate]; savefig(string(n.r.resultsDir,"mpc.",_pretty_defaults[:format])) end
  return h
end
"""
tp=tPlot(n,r,idx)
tp=tPlot(n,r,idx,tp;(:append=>true))
# plot the optimization times
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/11/2017, Last Modified: 6/28/2017 \n
--------------------------------------------------------------------------------------\n
"""
function tPlot(n,idx::Int64,args...;kwargs...);
  r = n.r.ocp

  kw = Dict(kwargs);
  # check to see if user would like to add to an existing plot
  if !haskey(kw,:append); append=false;
  else; append = get(kw,:append,0);
  end
  if !append; tp=plot(0,leg=:false); else tp=args[1]; end

  if idx > length(r.dfsOpt[:tSolve])
      warn("Cannot plot idx = ", idx, " because length(r.dfsOpt[:tSolve]) = ", length(r.dfsOpt[:tSolve]), ". \n
            reducing idx in tPlot!().")
      idx = length(r.dfsOpt[:tSolve])
  end

  # check to see if user would like to label legend
  if !haskey(kw,:legend);legend_string="";
  else; legend_string = get(kw,:legend,0);
  end

  # to avoid a bunch of jumping around in the simulation
	idx_max = length(r.dfsOpt[:tSolve])
	if (idx_max<10); idx_max=10 end

  if idx > 1
    scatter!(1:idx-1,r.dfsOpt[:tSolve][1:idx-1],marker=_pretty_defaults[:opt_marker][3],label=string(legend_string,"Previous Times"))
  end
  if isequal(r.dfsOpt[:status][idx],:Optimal)
    scatter!((idx,r.dfsOpt[:tSolve][idx]),marker=_pretty_defaults[:opt_marker][1],label=string(legend_string,"Optimal"))
  else
    scatter!((idx,r.dfsOpt[:tSolve][idx]),marker=_pretty_defaults[:opt_marker][2],label=string(legend_string,string(r.dfsOpt[:status][idx])) )
  end

  plot!(1:length(r.dfsOpt[:tSolve]),n.mpc.v.tex*ones(length(r.dfsOpt[:tSolve])),line=_pretty_defaults[:limit_lines][2],leg=:true,label="real-time threshhold",leg=:topright)

	ylims!(0,max(n.mpc.v.tex*1.2, maximum(r.dfsOpt[:tSolve])))
  xlims!(1,length(r.dfsOpt[:tSolve]))
	yaxis!("Optimization Time (s)")
	xaxis!("Evaluation Number")
  plot!(size=_pretty_defaults[:size]);
	if !_pretty_defaults[:simulate] savefig(string(n.r.resultsDir,"tplot.",_pretty_defaults[:format])) end

	return tp
end

"""
--------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 7/04/2017, Last Modified: 7/04/2017 \n
--------------------------------------------------------------------------------------\n
"""
function optPlot(n)
  L = length(n.r.ocp.dfsOpt[:tSolve])
  temp = [n.r.ocp.dfsOpt[:objVal][jj] for jj in 1:L]
  val = [idx for tempM in temp for idx=tempM]
  opt = plot(1:L,val)
  yaxis!("Objective Function Values")
	xaxis!("Evaluation Number")
  savefig(string(n.r.resultsDir,"optPlot.",_pretty_defaults[:format]))
  return opt
end
