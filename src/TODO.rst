TODO

==================
# High Priority #
==================
1) check if variable are getting overwritten!
2) configure calls define solver and that calls ocp def, so when you call define solver after configure this might be and issue, may just make it slow because you are calling ocpdef twice
3) why did it slow down!!? implicit?
4) start benchmarking so we can see the speed changes on a graph
5) do not export defineSolver!
===================
# Medium Priority #
===================
18) add ! to allPlots()
20) when integrating using :tm :trapezoidal method and minimizing the square of something it may oscillate
22) provide an initial guess
25) parse_DiffEq works OK except when in module
26) what does n.s.reset do
27) show difference in solution for explicit and implicit in Brachistochrone

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
