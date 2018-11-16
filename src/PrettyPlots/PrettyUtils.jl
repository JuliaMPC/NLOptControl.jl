#lc=distinguishable_colors(10);
const _pretty_defaults = Dict(
    :size         =>(700,700),   # size of figure
    :format       =>"png",       # format for output plots
    :plant        =>false,       # bool for simulation of plant as well as MPC
    :X0p          =>false,       # bool for simulation of X0p
    :X0p_marker  =>(:circle,:red,5.0,1.0),
    :plantOnly    =>false,       # bool for plotting only the plant response
    :polyPts      =>true,        # bool for plotting the individual polynomials
    :simulate     =>false,       # bool for simulation
    :L            =>100,         # format for output plots
    :mpc_lines    =>[(2.7,:blue,:solid),(2.7,:red,:solid),(2.7,:darkgreen,:solid),
                     (2.7,:yellow,:solid),(2.7,:black,:solid),(2.7,:firebrick,:solid),
                     (2.7,:gold,:solid),(2.7,:greenyellow,:solid),(2.7,:hotpink,:solid),
                     (2.7,:steelblue,:solid),(2.7,:aqua,:solid),(2.7,:teal,:solid),
                     (2.7,:blue,:dash),(2.7,:red,:dash),(2.7,:darkgreen,:dash),
                     (2.7,:yellow,:dash),(2.7,:black,:dash),(2.7,:firebrick,:dash),
                     (2.7,:gold,:dash),(2.7,:greenyellow,:dash),(2.7,:hotpink,:dash),
                     (2.7,:steelblue,:dash),(2.7,:aqua,:dash),(2.7,:teal,:dash)],
    :mpc_markers  =>(:circle,:blueviolet,0.0,1.0),
    :plant_lines  =>[(4.3,:darkgreen,:solid),(4.3,:darkgreen,:dash),(4.3,:black,:solid),(4.3,:black,:dash)],
    :limit_lines  =>[(2.0,:deeppink1,:solid),(2.0,:deeppink1,:dash),(2.0,:blueviolet,:solid),(2.0,:blueviolet,:dash)],
    :opt_marker  => [(:circle,:green,3.0),(:rect,:red,3.0),(:circle,:black,3.0)],

    # for VehicleModel.jl
    :vehicle_marker  =>(:black,:rect,5.0,1.0),
    :vehicle_fill    =>(0,1,:black),
    :goal_marker     =>[(:circle,:springgreen,10.0,1.0),(:circle,:green,10.0,1.0)],
    :goal_fill       =>[(0,1,:springgreen),(0,1,:green)],
    :goal_line       =>(1.0,0.0,:solid,:green),
    :obstacle_marker =>(:circle,:red,10.0,1.0),
    :obstacle_fill   =>(0,1.0,:red),
    :obstacle_line   =>(1.0,0.0,:solid,:red),
    :tire_force_lims =>(0,12000),
    :ax_lims         =>(-5,2),
    :vehicle_width   =>1.9,
    :vehicle_length  =>3.3,
    :lidar_line      =>(1.0,0.0,:solid,:red),
    :lidar_fill      =>(0,0.2,:yellow),
)

"""
plotSettings(;(:simulate=>true))
plotSettings(;(:mpc_lines   =>[(4.0,:blueviolet,:solid),(4.0,:blueviolet,:dash)]),
              (:plant_lines =>[(3.0,:darkgreen,:solid),(3.0,:darkgreen,:solid)]),
              (:limit_lines =>[(2.0,:turquoise,:solid),(2.0,:violet,:solid),(2.0,:orchid,:solid),(2.0,:darksalmon,:solid)]),
              (:size=>(1000,1000))
             )
plotSettings(;(:obstacle_marker =>(:circle,:red,10.0,1.0)),
                  (:goal_marker=>(:circle,:green,10.0,1.0))
             )
 plotSettings(;(:mpc_lines =>[(4.0,:blueviolet,:solid)]),(:size=>(700,700)))
# this function can be called to modify the default settings in PrettyPlots.jl
# check out Colors.jl for an incredible amount of colors!
"""
function plotSettings(;kwargs...)
  kw = Dict(kwargs);
  for (key,value) in kw
    if haskey(_pretty_defaults,key)
      _pretty_defaults[key]=value
    else
      error("Unknown key: ", kw)
    end
  end
end
"""
currentSettings()
# show the current plot settings
"""
function currentSettings()
  _pretty_defaults
end
