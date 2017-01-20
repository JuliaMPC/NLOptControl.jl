using DataFrames
@unpack numStates, lengthControlVector = nlp
@unpack Ni = ps
x0_con = c[1];
xf_con = c[2];
dyn_con = c[3];

x0c = DataFrame(state = 1:numStates, x0_con = getdual(x0_con));
xfc = DataFrame(state = 1:numStates, xf_con = getdual(xf_con));
bc = join(x0c,xfc, on=:step)

xc=Matrix{DataFrame}(Ni,numStates) # shows each interval and each state
for int in 1:Ni
  for st in 1:numStates
    xc[int,st] = DataFrame(i = 1:lengthControlVector, dyn_con = getdual(dyn_con)[int][:,st]);
  end
end
