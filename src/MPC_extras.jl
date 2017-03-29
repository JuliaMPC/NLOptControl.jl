# note this cannot be defined in the module because it uses a fucntion that is defined afterwards
"""
autonomousControl(mdl,n,r,s,params);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/1/2017, Last Modified: 3/11/2017 \n
--------------------------------------------------------------------------------------\n
"""
function autonomousControl(mdl,n,r,s,params)
 for st in 1:n.numStates   # update states based off of n.mpc.X0p
   if any(!isnan(n.X0_tol[st]))
     JuMP.setRHS(r.x0_con[st,1], (n.mpc.X0p[st]+n.X0_tol[st]));
     JuMP.setRHS(r.x0_con[st,2],-(n.mpc.X0p[st]-n.X0_tol[st]));
  else
    JuMP.setRHS(r.x0_con[st],n.mpc.X0p[st]);
   end
 end
 optimize(mdl,n,r,s)
end
