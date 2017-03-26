using NLOptControl, JuMP, Parameters, VehicleModels
main_dir=pwd();

# initialize
n = NLOpt(); s = Settings();
# define
pa = VparaKB(x0_=0.);  @unpack_VparaKB pa # vehicle parameters
X0 = [x0_,y0_,psi0_,u0_];
XF = [NaN,NaN,NaN,NaN];
XL = [x_min,y_min,psi_min,u_min];
XU = [x_max,y_max,psi_max,u_max];
CL = [sa_min,ax_min];
CU = [sa_max,ax_max];
n = define(n,stateEquations=KinematicBicycle,numStates=4,numControls=2,X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU)

# build
n = configure(n,Ni=2,Nck=[15,10];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 4.0))
defineSolver(n,solver=:IPOPT)
mdl = build(n);

# addtional information
names = [:x,:y,:psi,:ux];
descriptions = ["X (m)","Y (m)","Yaw Angle (rad)","Longitudinal Velocity (m/s)"];
stateNames(n,names,descriptions)
names = [:sr,:jx];
descriptions = ["Steering Angle (rad)","Longitudinal Acceleration (m/s^2)"];
controlNames(n,names,descriptions);

mXL=Any[false,false,false,false];mXU=Any[false,false,false,-1];  # set to false if you don't want to taper that side
linearStateTolerances(n;mXL=mXL,mXU=mXU);

# setup OCP
params = [pa];   # vehicle parameters
n,r = OCPdef(mdl,n,s,params);
x_ref = 10; y_ref = 100; # define target
@NLobjective(mdl, Min, (r.x[end,1]-x_ref)^2 + (r.x[end,2]-y_ref)^2);

# solve
optimize(mdl,n,r,s)

# post process
using PrettyPlots, Plots
gr();
allPlots(n,r,s,1)
