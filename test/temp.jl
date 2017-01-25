using NLOptControl
using JuMP
using Ipopt
using Parameters
using Plots
gr()

steps = 20; tests = 5;
global t_solve = Array(Float64, steps, tests)
global abs_error = Array(Float64, steps, tests)
global status = Array(Symbol, steps, tests)
global N = [i+1 for i in 1:steps]
global labels = Array(Symbol,tests)
##################################
# Define NLOptControl problem
##################################
function BrysonDenham{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2}) # dynamic constraint equations
  if n.integrationMethod==:tm
    L = size(x)[1];
  else
    L = size(x)[1]-1;
  end
  dx = Array(Any,L,n.numStates)
  dx[:,1] =  @NLexpression(mdl, [j=1:L], x[j,2] );
  dx[:,2] =  @NLexpression(mdl, [j=1:L], u[j,1] );
  return dx
end
L=1/9;
for q in 1:tests
    for j in 1:steps
        for m in 1:2
            if m==1
                global n = NLOpt();
                n = define(n,stateEquations=BrysonDenham,numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,-Inf],XU=[L,Inf],CL=[-Inf],CU=[Inf])
                if q==1
                    n = configure(n::NLOpt,Ni=1,Nck=[N[j]];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 1.0))
                elseif q==2
                    n = configure(n::NLOpt,N=N[j];(:integrationMethod => :tm),(:integrationScheme => :bkwEuler),(:finalTimeDV => false),(:tf => 1.0))
                elseif q==3
                    n = configure(n::NLOpt,N=N[j];(:integrationMethod => :tm),(:integrationScheme => :trapezoidal),(:finalTimeDV => false),(:tf => 1.0))
                elseif q==4
                    n = configure(n::NLOpt,Ni=2,Nck=[N[j],N[j]];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 1.0))
                else
                    n = configure(n::NLOpt,Ni=4,Nck=[N[j],N[j],N[j],N[j]];(:integrationMethod => :ps),(:integrationScheme => :lgrExplicit),(:finalTimeDV => false),(:tf => 1.0))
                end
                ##################################
                # Define JuMP problem
                ##################################
                global mdl = Model(solver = IpoptSolver(print_level=1));
                n,x,u,c=OCPdef(mdl,n)
                global obj = integrate(mdl,n,u[:,1];C=0.5,(:variable=>:control),(:integrand=>:squared))
                @NLobjective(mdl, Min, obj)
                solve(mdl)
            else
                t1 = time();
                status[j,q] = solve(mdl);
                t2 = time();
                t_solve[j,q] = t2 - t1;
                abs_error[j,q]=abs(getvalue(obj) - 4/(9*L))/(4/(9*L))*100;
            end
        end
    end
    if n.integrationMethod==:tm
        labels[q] = string(n.integrationScheme);
    else
        labels[q] = string(n.integrationScheme, " with Ni = ",n.Ni);
    end
end
