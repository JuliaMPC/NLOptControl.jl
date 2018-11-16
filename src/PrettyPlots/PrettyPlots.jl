module PrettyPlots

using Colors
using DataFrames
using Plots
using CSV
#import Plots.gr
#import Plots.pyplot
import Plots.xlims!, Plots.ylims!, Plots.plot
gr(); # default backend
#pyplot()

include("../Base.jl")
using .Base

include("PrettyUtils.jl")
include("NLOptControl_plots.jl")
include("MPCPlots.jl")
#include("VehicleModels_plots.jl")

export
    # PrettyUtils.jl
    minDF,
    maxDF,
    plotSettings,
    _pretty_defaults,
    currentSettings,

    # NLOptControl.jl plots
    statePlot,
    controlPlot,
    costatesPlot,
    costatesPlots,
    allPlots,
    adjust_axis,

    # MPC plots
    mpcPlot,
    tPlot,
    optPlot,

    # Plots.jl exported functions
    xlims!,
    ylims!,
    plot
    #Plots.gui(),
    #Plots.pyplot(),
    #Plots.gr(),
    #Plots.pgfplots()

end # module
