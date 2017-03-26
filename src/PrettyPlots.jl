"""
resultsDir(results_dir)
# removes results folder and creates a new one
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function resultsDir(results_dir::String)
	if isdir(results_dir)
		rm(results_dir; recursive=true)
		print("\n The old results have all been deleted! \n \n")
	end
	mkdir(results_dir)
	cd(results_dir)
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 3/11/2017 \n
--------------------------------------------------------------------------------------\n
"""
function allPlots(n::NLOpt,r::Result,s::Settings,idx::Int64)
  stp = [statePlot(n,r,s,idx,st) for st in 1:n.numStates];
  ctp = [controlPlot(n,r,s,idx,ctr) for ctr in 1:n.numControls];
  all = [stp;ctp];
  h = plot(all...,size=(s.s1,s.s1));
  if !s.simulate; savefig(string(r.results_dir,"main.",s.format)) end
  return h
end

"""
stp=statePlot(n,r,s,idx,st);
stp=statePlot(n,r,s,idx,st;(:legend=>"test1"));
stp=statePlot(n,r,s,idx,st,stp;(:append=>true));
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 3/11/2017 \n
--------------------------------------------------------------------------------------\n
"""
function statePlot(n::NLOpt,r::Result,s::Settings,idx::Int64,st::Int64,args...;kwargs...)
  kw = Dict(kwargs);

  # check to se if user would like to add to an existing plot
  if !haskey(kw,:append); kw_ = Dict(:append => false); append = get(kw_,:append,0);
  else; append = get(kw,:append,0);
  end
  if !append; stp=plot(0,leg=:false); else stp=args[1]; end

  # check to see if user would like to plot limits
  if !haskey(kw,:lims); kw_ = Dict(:lims => true); lims = get(kw_,:lims,0);
  else; lims = get(kw,:lims,0);
  end

  # check to see if user would like to label legend
  if !haskey(kw,:legend); kw_ = Dict(:legend => ""); legend_string = get(kw_,:legend,0);
  else; legend_string = get(kw,:legend,0);
  end

	if r.dfs[idx]!=nothing
  	t_vec=linspace(0.0,max(5,round(r.dfs[end][:t][end]/5)*5),s.L);
	else
		t_vec=linspace(r.dfs_plant[1][:t][1],max(5,round(r.dfs_plant[end][:t][end]/5)*5),s.L);
	end

  # plot the limits
  if lims
		if n.linearStateTol[st]
			if !isinf(n.XL[st]);plot!(r.t_st,n.XL_var[st,:],w=s.lw1,label=string(legend_string,"min"));end
    	if !isinf(n.XU[st]);plot!(r.t_st,n.XU_var[st,:],w=s.lw1,label=string(legend_string,"max"));end
		else
    	if !isinf(n.XL[st]);plot!(t_vec,n.XL[st]*ones(s.L,1),w=s.lw1,label=string(legend_string,"min"));end
    	if !isinf(n.XU[st]);plot!(t_vec,n.XU[st]*ones(s.L,1),w=s.lw1,label=string(legend_string,"max"));end
		end
  end

  # plot the values
	if r.dfs[idx]!=nothing
  	plot!(r.dfs[idx][:t],r.dfs[idx][n.state.name[st]],w=s.lw1,label=string(legend_string,"mpc"));
	end
  if s.MPC
		# values
		temp = [r.dfs_plant[jj][n.state.name[st]] for jj in 1:idx];
	  vals=[idx for tempM in temp for idx=tempM];

		# time
		temp = [r.dfs_plant[jj][:t] for jj in 1:idx];
		time=[idx for tempM in temp for idx=tempM];

    plot!(time,vals,w=s.lw2,label=string(legend_string,"plant"));
  end
  adjust_axis(xlims(),ylims());
	xlims!(t_vec[1],t_vec[end]);
  plot!(size=(s.s1,s.s1));
  yaxis!(n.state.description[st]); xaxis!("time (s)");
  if !s.simulate; savefig(string(r.results_dir,n.state.name[st],".",s.format)); end
  return stp
end

"""
stp=statePlot(n,r,s,idx,st1,st2);
stp=statePlot(n,r,s,idx,st1,st2;(:legend=>"test1"));
stp=statePlot(n,r,s,idx,st1,st2,stp;(:append=>true),(:lims=>false));
# to compare two different states
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 3/11/2017 \n
--------------------------------------------------------------------------------------\n
"""
function statePlot(n::NLOpt,r::Result,s::Settings,idx::Int64,st1::Int64,st2::Int64,args...;kwargs...)
  kw = Dict(kwargs);

  # check to see if user would like to add to an existing plot
  if !haskey(kw,:append); kw_ = Dict(:append => false); append = get(kw_,:append,0);
  else; append = get(kw,:append,0);
  end
  if !append; stp=plot(0,leg=:false); else stp=args[1]; end

  # check to see if user would like to plot limits
  if !haskey(kw,:lims); kw_ = Dict(:lims => true); lims = get(kw_,:lims,0);
  else; lims = get(kw,:lims,0);
  end

  # check to see if user would like to label legend
  if !haskey(kw,:legend); kw_ = Dict(:legend => ""); legend_string = get(kw_,:legend,0);
  else; legend_string = get(kw,:legend,0);
  end

  # plot the limits
  # TODO check if all constraints are given
  if lims
    if !isinf(n.XL[st1]);plot!([n.XL[st1],n.XL[st1]],[n.XL[st2],n.XU[st2]],w=s.lw1,label=string(n.state.name[st1],"_min"));end
    if !isinf(n.XU[st1]);plot!([n.XU[st1],n.XU[st1]],[n.XL[st2],n.XU[st2]],w=s.lw1,label=string(n.state.name[st1],"_max"));end

    if !isinf(n.XL[st1]);plot!([n.XL[st1],n.XU[st1]],[n.XL[st2],n.XL[st2]],w=s.lw1,label=string(n.state.name[st2],"_min"));end
    if !isinf(n.XU[st1]);plot!([n.XL[st1],n.XU[st1]],[n.XU[st2],n.XU[st2]],w=s.lw1,label=string(n.state.name[st2],"_max"));end
  end

  # plot the values
	if r.dfs[idx]!=nothing
		plot!(r.dfs[idx][n.state.name[st1]],r.dfs[idx][n.state.name[st2]],w=s.lw1,label=string(legend_string,"mpc"));
	end

  if s.MPC
		# values
		temp = [r.dfs_plant[jj][n.state.name[st1]] for jj in 1:idx];
		vals1=[idx for tempM in temp for idx=tempM];

		# values
		temp = [r.dfs_plant[jj][n.state.name[st2]] for jj in 1:idx];
		vals2=[idx for tempM in temp for idx=tempM];

		plot!(vals1,vals2,w=s.lw2,label=string(legend_string,"plant"));
  end
  adjust_axis(xlims(),ylims());
  plot!(size=(s.s1,s.s1));
  xaxis!(n.state.description[st1]);
  yaxis!(n.state.description[st2]);
  if !s.simulate savefig(string(r.results_dir,n.state.name[st1],"_vs_",n.state.name[st2],".",s.format)); end
  return stp
end

"""
ctrp=controlPlot(n,r,s,idx,ctr);
ctrp=controlPlot(n,r,s,idx,ctr,ctrp;(:append=>true));
# to plot control signals
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 3/11/2017 \n
--------------------------------------------------------------------------------------\n
"""
function controlPlot(n::NLOpt,r::Result,s::Settings,idx::Int64,ctr::Int64,args...;kwargs...)
  kw = Dict(kwargs);

  # check to see if user would like to add to an existing plot
  if !haskey(kw,:append); kw_ = Dict(:append => false); append = get(kw_,:append,0);
  else; append = get(kw,:append,0);
  end
  if !append; ctrp=plot(0,leg=:false); else ctrp=args[1]; end

  # check to see if user would like to plot limits
  if !haskey(kw,:lims); kw_ = Dict(:lims => true); lims = get(kw_,:lims,0);
  else; lims = get(kw,:lims,0);
  end

  # check to see if user would like to label legend
  if !haskey(kw,:legend); kw_ = Dict(:legend => ""); legend_string = get(kw_,:legend,0);
  else; legend_string = get(kw,:legend,0);
  end

	if r.dfs[idx]!=nothing
		t_vec=linspace(0.0,max(5,round(r.dfs[end][:t][end]/5)*5),s.L);
	else
		t_vec=linspace(r.dfs_plant[1][:t][1],max(5,round(r.dfs_plant[end][:t][end]/5)*5),s.L);
	end

  # plot the limits
  if lims
    if !isinf(n.CL[ctr]); plot!(t_vec,n.CL[ctr]*ones(s.L,1),w=s.lw1,label="min"); end
    if !isinf(n.CU[ctr]);plot!(t_vec,n.CU[ctr]*ones(s.L,1),w=s.lw1,label="max");end
  end

  # plot the values
	if r.dfs[idx]!=nothing
  	plot!(r.dfs[idx][:t],r.dfs[idx][n.control.name[ctr]],w=s.lw1,label=string(legend_string,"mpc"));
	end
  if s.MPC
		# values
		temp = [r.dfs_plant[jj][n.state.name[st]] for jj in 1:idx];
	  vals=[idx for tempM in temp for idx=tempM];

		# time
		temp = [r.dfs_plant[jj][:t] for jj in 1:idx];
		time=[idx for tempM in temp for idx=tempM];

		plot!(time,vals,w=s.lw2,label=string(legend_string,"plant"));
  end
  adjust_axis(xlims(),ylims());
	xlims!(t_vec[1],t_vec[end]);
  plot!(size=(s.s1,s.s1));
  yaxis!(n.control.description[ctr]);	xaxis!("time (s)");
	if !s.simulate savefig(string(r.results_dir,n.control.name[ctr],".",s.format)) end
  return ctrp
end

function adjust_axis(x_lim,y_lim)

	# scaling factors
	al_x = [0.05, 0.05]; # x axis (low, high)
	al_y = [0.05, 0.05];  # y axis (low, high)

	# additional axis movement
	if x_lim[1]==0.0; a=-1; else a=0; end
	if x_lim[2]==0.0; b=1;  else b=0; end
	if y_lim[1]==0.0; c=-0.01; else c=0; end
	if y_lim[2]==0.0; d=1;  else d=0; end

	xlim = Float64[0,0]; ylim = Float64[0,0];
	xlim[1] = x_lim[1]+x_lim[1]*al_x[1]+a;
	xlim[2] = x_lim[2]+x_lim[2]*al_x[2]+b;
	ylim[1] = y_lim[1]+y_lim[1]*al_y[1]+c;
	ylim[2] = y_lim[2]+y_lim[2]*al_y[2]+d;

	xlims!((xlim[1],xlim[2]))
	ylims!((ylim[1],ylim[2]))
end
