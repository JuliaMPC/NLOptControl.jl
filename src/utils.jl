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
