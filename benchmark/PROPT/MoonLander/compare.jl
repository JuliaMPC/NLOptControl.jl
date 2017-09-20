using DataFrames
using Plots
pgfplots()


# load optimal solution



# import NLOptControl.jl results

# import PROPT results



# compare PROPT and NLOptControl.jl


# plot the results
results_dir=string("benchmark1_",c.m.name,"/")
resultsDir!(n;results_name=results_dir);

# optimization times
optP=scatter(1:10,t_optJ,label="julia")
scatter!(1:10,t_optM,label="MATLAB")
yaxis!("Optimization Time (s)");
xaxis!("Evaluation Number");
title!(string("MATLAB ave. = ",round(aveM,2) ,"julia ave. = ",round(aveJ,2) ))
savefig(string(n.r.results_dir,"/","optPlot",".",:svg));


scatter(t_opt,h_opt,label="h_opt")
scatter!(t_opt,v_opt,label="v_opt")
scatter!(t_opt,u_opt,label="u_opt")
