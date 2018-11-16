using DataFrames
using Plots
using ImplicitEquations
using VehicleModels

# Set The Problem Up
main_dir="/home/febbo/Documents/workspace/OCP"
results_name = "testing_clp76"; # Name Folder to Store Comparison Results in
results_dir = string(main_dir,"/results/",results_name)
warn("you may have the wrong senario!")
senario = "s1.jl"  #TODO save this with the actual test data!

# define vehicle parameters; has constraints on vehicle model for optimization
pa = Vpara();   # initialize parameter set
@unpack_Vpara pa # unpack parameters

include(string(main_dir,"/initialize.jl"))


# Choose Results to Compare
k = 4 # number of results to compare
path = ["","","",""] # should match k
path[1] = "/home/febbo/Documents/workspace/OCP/results/testing_clp76eulerIPOPT";
path[2] = "/home/febbo/Documents/workspace/OCP/results/testing_clp76trap2IPOPT";
path[3] = "/home/febbo/Documents/workspace/OCP/results/testing_clp76eulerKNITRO";
path[4] = "/home/febbo/Documents/workspace/OCP/results/testing_clp76KNITRO";

# Plot Settings
label_string = ["Euler & IPOPT", "Trap. & IPOPT","Euler & KNITRO", "Trap. & KNITRO"];

if isdir(results_dir)
		rm(results_dir; recursive=true)
		print("\n The old results have all been deleted! \n \n")
end

dfs=data_sets(k,"STATES")  # put the data from each test into a DataFrame
dfs_opt=data_sets(k,"OPT_info")
mkdir(results_dir)
cd(results_dir)

sa=saplot(dfs,label_string,false,pa)
sr=srplot(dfs,label_string,false,pa)
yaw=yawplot(dfs,label_string,false,pa)
yawr=yawrplot(dfs,label_string,false,pa)

longv=longvplot(dfs,label_string,false,pa)
ax=axplot(dfs,label_string,false,pa)
jx=jxplot(dfs,label_string,false,pa)

latv=latvplot(dfs,label_string,false,pa)

pp=pplot(dfs,1,label_string,false,false,obs_data,s_data,pa)

lt=ltplot(dfs,label_string,false,pa)
vt=vtplot(dfs,label_string,false,pa)

# make a figure with all states and controls
plot(longv, ax, jx, sa, sr, yaw, yawr, latv, pp, size=(1800,1200))
savefig("main.png")

# make a figure with the tire forces
plot(lt, vt, size=(1800,1200),layout = grid(2,1) )
savefig("tire_f.png")

# make a simulation with all vehicles
#panim_fun(label_string)

optplot(dfs_opt,label_string,false)

tplot(dfs_opt,dfs,label_string,false)

cd(main_dir)
