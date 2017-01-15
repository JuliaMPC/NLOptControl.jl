macro OCPdef(mdl,nlp,ps)
#TODO figure out macros!
show(args)
mdl =esc(mdl)
show(mdl)
nlp = esc(nlp)
ps = ps
  #mdl = esc(mdl); nlp = esc(nlp); ps = esc(ps);
  #ss = [Symbol("xy$(st)") for st in 1:numStates];
  #for st in 1:numStates TODO figure out how to generate variable names
    #@variable(mdl, XL[st] <= xx1[1:sum(numStatePoints)] <= XU[st])
    #@variable(mdl, XL[st] <= $string(ss[st])[1:sum(numStatePoints)] <= XU[st])
    #anonymous2 =  @variable(mdl, [pts in 1:5])
  #end
print("\n",nlp,"\n" )
print("\n", mdl,"\n" )
  # inequality constraints and design variable definitions
  @unpack numStates, numStatePoints,  XL, XU = nlp
  @variable(mdl, x[1:sum(numStatePoints),1:numStates]);
  for st in 1:numStates
    for j in 1:sum(numStatePoints)
      setlowerbound(x[j,st], XL[st])
      setupperbound(x[j,st], XU[st])
    end
  end

  @unpack numControls, numControlPoints, CL, CU = nlp
  @variable(mdl, u[1:sum(numControlPoints),1:numControls]);
  for ctr in 1:numControls
    for j in 1:sum(numControlPoints)
      setlowerbound(u[j,ctr], CL[ctr])
      setupperbound(u[j,ctr], CU[ctr])
    end
  end



end
