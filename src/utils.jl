"""
# integrating JuMP variables
Expr = integrate(mdl,n,u;(:mode=>:control))
Expr = integrate(mdl,n,u,idx=1;C=0.5,(:variable=>:control),(:integrand=>:squared))
Expr = integrate(mdl,n,r.u[:,1];D=1.2,(:variable=>:control),(:integrand=>:squared),(:integrandAlgebra=>:subtract))

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 1/27/2017 \n
--------------------------------------------------------------------------------------\n
"""
function integrate(mdl::JuMP.Model,n::NLOpt,V::Array{JuMP.Variable,1}, args...; C::Float64=1.0,D=0.0,kwargs...)
  kw = Dict(kwargs);
  if !haskey(kw,:integrand); kw_ = Dict(:integrand => :default); integrand = get(kw_,:integrand,0);
  else; integrand = get(kw,:integrand,0);
  end
  if !haskey(kw,:integrandAlgebra); kw_ = Dict(:integrandAlgebra => :default); integrandAlgebra = get(kw_,:integrandAlgebra,0);
  else; integrandAlgebra = get(kw,:integrandAlgebra,0);
  end

  if integrandAlgebra ==:subtract
    if n.integrationMethod==:tm
      if integrand == :default      # integrate V
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl,C*sum((V[j+1]-D)*n.tf/(n.N) for j = 1:n.N));  #TODO fix this.. there is an extra dv here for control, but it does not effect solution
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl,C*sum(0.5*((V[j]-D+V[j+1]-D)*n.tf/(n.N) for j = 1:n.N)));
        end
      elseif integrand == :squared # integrate V^2
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl, C*sum(((V[j+1]-D)^2)*n.tf/(n.N) for j = 1:n.N));
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl, C*sum(0.5*((V[j]-D)^2+(V[j+1]-D)^2)*n.tf/(n.N) for j = 1:n.N));
        end
      else
        error("\n Check :integrand \n")
      end
    elseif n.integrationMethod==:ps
      if !haskey(kw,:mode); kw_ = Dict(:mode => :quadrature); mode = get(kw_,:mode,0);
      else; mode  = get(kw,:mode,0);
      end
      variable = get(kw,:variable,0);
      if variable == :state; Nck_cum  = [0;cumsum(n.Nck+1)];
      elseif variable == :control; Nck_cum = [0;cumsum(n.Nck)];
      else; error("\n Set the variable to either (:variable => :state) or (:variable => :control). \n")
      end

      if mode == :quadrature  #TODO recalculate ws based off of time
        if integrand == :default      # integrate V
          @NLexpression(mdl, temp[int=1:n.Ni], ((n.tf-n.t0)/2)*sum((n.ωₛ[int])[j]*((V[Nck_cum[int]+1:Nck_cum[int+1]])[j]-D) for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, C*sum(temp[int] for int = 1:n.Ni));
        elseif integrand == :squared # integrate V^2
          @NLexpression(mdl, temp[int=1:n.Ni],((n.tf-n.t0)/2)*C*sum((n.ωₛ[int])[j]*((V[Nck_cum[int]+1:Nck_cum[int+1]])[j]-D)*((V[Nck_cum[int] + 1:Nck_cum[int+1]])[j]-D) for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, sum(temp[int] for int = 1:n.Ni));
        else
          error("\n Check :integrand \n")
        end
      elseif mode == :LGRIM# TODO add in option to allow for integration using IMatrix
          error("\n Not implemented yet!! \n")
      end
    end
  elseif integrandAlgebra==:default
    if n.integrationMethod==:tm
      if integrand == :default      # integrate V
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl,C*sum(V[j+1]*n.tf/(n.N) for j = 1:n.N));  #TODO fix this.. there is an extra dv here for control, but it does not effect solution
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl,C*sum(0.5*(V[j]+V[j+1])*n.tf/(n.N) for j = 1:n.N));
        end
      elseif integrand == :squared # integrate V^2
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl, C*sum((V[j+1]^2)*n.tf/(n.N) for j = 1:n.N));
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl, C*sum(0.5*(V[j]^2+V[j+1]^2)*n.tf/(n.N) for j = 1:n.N));
        end
      else
        error("\n Check :integrand \n")
      end
    elseif n.integrationMethod==:ps
      if !haskey(kw,:mode); kw_ = Dict(:mode => :quadrature); mode = get(kw_,:mode,0);
      else; mode  = get(kw,:mode,0);
      end
      variable = get(kw,:variable,0);
      if variable == :state; Nck_cum  = [0;cumsum(n.Nck+1)];
      elseif variable == :control; Nck_cum = [0;cumsum(n.Nck)];
      else; error("\n Set the variable to either (:variable => :state) or (:variable => :control). \n")
      end

      if mode == :quadrature  #TODO recalculate ws based off of time
        if integrand == :default      # integrate V
          @NLexpression(mdl, temp[int=1:n.Ni], ((n.tf-n.t0)/2)*sum((n.ωₛ[int])[j]*(V[Nck_cum[int]+1:Nck_cum[int+1]])[j] for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, C*sum(temp[int] for int = 1:n.Ni));
        elseif integrand == :squared # integrate V^2
          @NLexpression(mdl, temp[int=1:n.Ni],((n.tf-n.t0)/2)*C*sum((n.ωₛ[int])[j]*(V[Nck_cum[int]+1:Nck_cum[int+1]])[j]*(V[Nck_cum[int] + 1:Nck_cum[int+1]])[j] for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, sum(temp[int] for int = 1:n.Ni));
        else
          error("\n Check :integrand \n")
        end
      elseif mode == :LGRIM# TODO add in option to allow for integration using IMatrix
          error("\n Not implemented yet!! \n")
      end
    end
  end
  return Expr
end

"""
mdl=build(n)
# build JuMP model
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 2/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function build(n::NLOpt) #TODO allow user to pass solver options
  if n.solver==:IPOPT
    mdl=Model(solver=IpoptSolver());
  else n.solver==:KNITRO
    mdl=Model(solver=KnitroSolver());
  end
  return mdl
end

"""
optimize(mdl,n,r)
# solves JuMP model and saves optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/6/2017, Last Modified: 2/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function optimize(mdl::JuMP.Model, n::NLOpt, r::Result)
  t1 = time(); status = JuMP.solve(mdl); t2 = time();
  push!(r.status,status);
  push!(r.t_solve,(t2 - t1));
  push!(r.obj_val, getobjectivevalue(mdl));
  r.eval_num=length(r.status);
  postProcess(n,r);
end


"""
# funtionality to save constraint data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 2/8/2017 \n
--------------------------------------------------------------------------------------\n
"""
type Constraint
  name::Vector{Any}
  handle::Vector{Any}
  value::Vector{Any}
end

function Constraint()
  Constraint([],
             [],
             []);
end

function newConstraint(r::Result,handle,name::Symbol)
  initConstraint(r)
  constraint::Constraint = r.constraint
  push!(r.constraint.handle,handle)
  push!(r.constraint.name,name)
end

function initConstraint(r::Result)
  if r.constraint == nothing
    r.constraint = Constraint()
  end
end

function evalConstraints(n::NLOpt, r::Result)
  r.constraint.value=[];   # reset values
  for i = 1:length(r.constraint.handle)
    if r.constraint.name[i]==:dyn_con  # state constraits
      dfs=Vector{DataFrame}(n.numStates);
      con=DataFrame(step=1);
      for st in 1:n.numStates
        if n.integrationMethod==:ps
          temp=[getdual(r.constraint.handle[i][int][:,st]) for int in 1:n.Ni];
          vals=[idx for tempM in temp for idx=tempM];
          dfs[st]=DataFrame(step=1:sum(n.Nck);Dict(n.state.name[st] => vals)...);
        else
          dfs[st]=DataFrame(step=1:n.N;Dict(n.state.name[st] => getdual(r.constraint.handle[i][:,st]))...);
        end
        if st==1; con=dfs[st]; else; con=join(con,dfs[st],on=:step); end
      end
    else
      S=JuMP.size(r.constraint.handle[i])
      if length(S)==1
        con = DataFrame(step=1:length(r.constraint.handle[i]);Dict(r.constraint.name[i] => getdual(r.constraint.handle[i][:]))...);
      elseif length(S)==2
        dfs=Vector{DataFrame}(S[1]);
        con=DataFrame(step=1);
        for idx in 1:S[1]
          dfs[idx] = DataFrame(step=1:S[2];Dict(r.constraint.name[i] => getdual(r.constraint.handle[i][idx,:]))...);
          if idx==1; con=dfs[idx]; else; con=join(con,dfs[idx],on=:step); end
        end
      end
    end
    push!(r.constraint.value,con)
  end
end

"""
# funtionality for state data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 2/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
type State #TODO correlate these with JuMP variables
  name::Vector{Any}
  description::Vector{Any}
end

function State()
  State([],
        []);
end

function initStateNames(n::NLOpt)
  State([Symbol("x$i") for i in 1:n.numStates],
        [String("x$i") for i in 1:n.numStates]);
end

function stateNames(n::NLOpt,names,descriptions)
  n.state::State = State() # reset
  if length(names)!=n.numStates || length(descriptions)!=n.numStates
    error("\n Check sizes of names and descriptions \n")
  end
  for i in 1:n.numStates
    push!(n.state.name,names[i])
    push!(n.state.description,descriptions[i])
  end
end

"""
# funtionality for control data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 2/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
type Control #TODO correlate these with JuMP variables
  name::Vector{Any}
  description::Vector{Any}
end

function Control()
  Control([],
        []);
end

function initControlNames(n::NLOpt)
  Control([Symbol("u$i") for i in 1:n.numControls],
          [String("u$i") for i in 1:n.numControls]);
end

function controlNames(n::NLOpt,names,descriptions)
  n.control::Control = Control() # reset
  if length(names)!=n.numControls || length(descriptions)!=n.numControls
    error("\n Check sizes of names and descriptions \n")
  end
  for i in 1:n.numControls
    push!(n.control.name,names[i])
    push!(n.control.description,descriptions[i])
  end
end
