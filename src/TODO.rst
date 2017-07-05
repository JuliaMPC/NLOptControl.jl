TODO

==================
# High Priority #
==================
1) check if variable are getting overwritten!
2) configure calls define solver and that calls ocp def, so when you call define solver after configure this might be and issue, may just make it slow because you are calling ocpdef twice
3) why did it slow down!!? implicit?
4) start benchmarking so we can see the speed changes on a graph
5) look at the time horizon for the MAVs problem!!
6) add an option for tf max or maybe REDUCE IT AT LEAST!!!!
7) make a better test for MAVs
8) using NLExpr are we nesting extra expressions? maybe that is why it is slow?
9) upper acceleration constraint not met...
10) maybe there needs to be deepcopy() of NLexpr variables?
11) maybe it has to do with how the NLP is organized? we moved a bunch of constraints to a new part of the problem...
try making an additional function where constraints can be added then putting it into configure so that all of the constraints are added up
12) try to return dx[:,1] =
13) it is strange that when n.evalConstraints==true the solution can change!
14) put functions back in for vehicle model and compare the results
 dynamics do not need to be in expressions, just the tire model for the obj function
15) NaN or -Inf and Inf?
===================
# Medium Priority #
===================
18) add ! to allPlots()
20) when integrating using :tm :trapezoidal method and minimizing the square of something it may oscillate
22) provide an initial guess
25) parse_DiffEq works OK except when in module
26) what does n.s.reset do
27) show difference in solution for explicit and implicit in Brachistochrone
28) WARNING: Dual solution not available. Check that the model was properly solved and no integer variables are present.
don;t try to extract if the solution is not feasible
29) add a test for constraints!()
30) why does the problem solve so poorly using setsolver?

=================
# Low Priority #
=================
1) warn user if final time constraint is active
6) try other linear solvers with Ipopt
14) make functionality to easily compare all integration schemes
15) implement method with variable time steps for :tm methods, that is how it was before!
try to register functions with JuMP
17) change warning from Pkg.add(`` ...) to ""
19) fix pgfplots in mpcdocs
20) get a simple MPC example working?
21) get ride of the index: de=[:(v[j]),:(T[j]-1.625)]
22) make it stop recompiling VehicleModels when it runs PrettyPlots
23) see the number of constraints can be reduced for :ps methods
25) update notebook
30) if user only passes n and Expr to NLExpr automatically add n.r.x etc.
31) make this a oneliner
 FZ_rl_con=@NLconstraint(n.mdl, [j=1:n.numStatePoints], 0 <= FZ_RL[j] - Fz_min)
 newConstraint!(n,FZ_rl_con,:FZ_rl_con);
33) To debug KNITRO turn up the output level
34) Try to tune KNITRO
35) put an example of not reducing tf
36) NLcon and NLexpr have very similar code, consider merging
37) add in   newConstraint!(n,FZ_rl_con,:FZ_rl_con);
38) think about how to make the n object imutable unless it is modified in NLOPtCOntrol using one of it's functions
39) consider elliminating the NLexpressions and directly using NLconstraints for dynamics
40) eventually consider adding parameters so that obstacle avoidance constraints can be added in sequence
41)  @variable( mdl, 0.00001 <= dt[1:n.N] <= 0.2) #TODO allow for an varaible array of dts
41) put this warning in newConstraint
         error("\n For now, the constraints cannot be in this form: \n
        con=@NLconstraint(mdl,n.r.u[1,1]==param); \n
        Write it in array form: \n
          con=@NLconstraint(mdl,[i=1],n.r.u[i,1]==param); \n")
