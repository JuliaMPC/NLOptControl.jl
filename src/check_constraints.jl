using DataFrames

x0_con = c[1];
xf_con = c[2];
dyn_con = c[3];

x0c = DataFrame(state = 1:n.numStates, x0_con = getdual(x0_con));
xfc = DataFrame(state = 1:n.numStates, xf_con = getdual(xf_con));
bc = join(x0c,xfc, on=:state)

xc=Matrix{DataFrame}(n.Ni,n.numStates) # shows each interval and each state
for int in 1:n.Ni
  for st in 1:n.numStates
    xc[int,st] = DataFrame(i = 1:n.numPoints[int], dyn_con = getdual(dyn_con)[int][:,st]);
  end
end
