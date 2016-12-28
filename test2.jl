using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature
pyplot()
#default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

"""
# variables
* int = index for interval
* inx = index for collocation point
"""
# define problem
t0 = Float64(0); tf = Float64(10);  # TODO change and y to x
t = Array(linspace(t0,tf,100));
α₁ =  -0.3; α₂ = 3; α₃ = -8; α₄ =  7;
Nc = Int64(3); # number of collocation points in each interval
Ni = Int64(2);  # number of intervals
####################################
# perform analytical calcualtions
####################################
γ = Poly([α₁,α₂,α₁]); #TODO check on that imported binding warning
y = polyval(γ,t);

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,t[end]) - polyval(∫γ,t[1]);
C = Y - polyval(∫γ,t[end]); # constant of integration
∫y = polyval(∫γ,t) + C;

# evaluate the derivative
dγ = polyder(γ);
dy = polyval(dγ,t);

####################################
# construct polynomial approximation
####################################
τ, ω = gaussradau(Nc); # number of collocation points per interval

# break the problem up into multiple intervals
di, tm, t_data, ωₛ = create_intervals(t0,tf,Ni,Nc,τ,ω); # TODO probably get ride of some of the output

temp = zeros(Float64,1,Ni)
for int in 1:Ni
  temp[int] = di*int;
end
t_dataE = [t_data;temp];
# approximate state (or polynomial)
state_data = zeros(Float64,Nc+1,Ni);
for int in 1:Ni
  state_data[:,int] = polyval(γ,t_dataE[:,int]);
end
P = zeros(Float64,Nc+1,Ni);
for int in 1:Ni
  P[:,int] = interpolate_lagrange(t_dataE[:,int],t_dataE[:,int],state_data[:,int],Nc)
end

# approximate the integral
ζ = zeros(Float64,Nc,Ni); fτ = zeros(Float64,Nc,Ni);approx_int = Float64(0);
for idx in 1:Ni
  fτ[:,idx] = polyval(γ,t_data[:,idx]);
  ζ[:,idx] =  cumsum(ωₛ[:,idx].*fτ[:,idx],1)
  approx_int = approx_int + ζ[end,idx];
end
# approximate the derivative --> needed in defect constraints
dζ = zeros(Float64,Nc+1,Ni);fτE = zeros(Float64,Nc+1,Ni); D = zeros(Float64,Nc+1,Nc+1,Ni);
for int in 1:Ni
  fτE[:,int] = polyval(γ,t_dataE[:,int]);
  D[:,:,int] = poldif(t_dataE[:,int], 1)
  dζ[:,int] = D[:,:,int]*fτE[:,int]
end

#################
# post processing
#################
tF = zeros(Float64,Ni); yF =  zeros(Float64,Ni);
fp=plot(0,leg=:false);
plot!(t,y,label=string(γ),w=6)
for int in 1:Ni
  scatter!(t_dataE[1:end-1,int],P[1:end-1,int],markersize =10,markershape = :rect,label=string("collocation points for mesh interval ",int))
  tF[int] = t_dataE[end,int];
  yF[int] = P[end,int];
end
scatter!(tF,yF,markersize = 10,markershape = :star8,label=string("end points"))
xlims!(t0,tf*1.1)

dp=plot(0,leg=:false)
plot!(t,dy,label="derivative",w=6)
for int in 1:Ni
  scatter!(t_dataE[1:end-1,int],dζ[1:end-1,int],markersize =10,markershape = :rect,label=string("approximate derivative ",int))
  tF[int] = t_dataE[end,int];
  yF[int] = dζ[end,int];
end
scatter!(tF,yF,markersize = 10,markershape = :star8,label=string("end points"))
xlims!(t0,tf*1.1)

ip=plot(0,leg=:false)
plot!(t,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=6)
plot!(t_data,ζ,label=@sprintf("approx. integral = %0.3f",approx_int),line=(4,:dash))

plot(fp,ip,dp,layout=(3,1))


#=

temp = zeros(Float64,1,Ni)
for idx in 1:Ni
  temp[idx] = di*idx;
end
t_dataE = [t_data;temp];
# approximate state (or polynomial)
state_data = zeros(Float64,Nc+1,Ni);
for j in 1:Ni
  state_data[:,j] = polyval(γ,t_dataE[:,j]);
end

# Calculate the Differention Matrixes
L = zeros(Float64,Nc+1,Nc,Ni);
for j in 1:Ni
  L[:,:,j] = lagrange_basis_poly(t_dataE[1:end-1,j],t_dataE[:,j],Nc)  # TODO check -->may be wrong size [N]X[N+1]
end

P = zeros(Float64,Nc+1,Nc,Ni)
for j in 1:Ni
  for idx in 1:Nc+1
      if idx==1
          P[idx,:,j] = state_data[idx,j]*L[idx,:,j]
      else
          P[idx,:,j] = state_data[idx,j]*L[idx,:,j] + P[idx-1,:,j]
      end
  end
end

plot(t,y,label=string(γ),w=6)
plot!(t_dataE[1:end-1,1],P[1,:,1],w=8,label="P0")
plot!(t_dataE[1:end-1,1],P[2,:,1],w=8,label="P1")
plot!(t_dataE[1:end-1,1],P[3,:,1],w=8,label="P2")
plot!(t_dataE[1:end-1,1],P[4,:,1],line=(8,:dash,:white),label="P3",background_color=RGB(0.2,0.2,0.2))
scatter!(t_data[:,1],state_data[:,1],markersize =12,markershape = :hexagon,label="data")
=#
#ylims!(0,10)
