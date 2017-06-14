TODO

# today
1) use the DiffEq macro
2) modify examples
  - get ride of using Plots etc.
3) start a julia notebook
4) get juliaBox working

# tomorrow
1) finish julia notebook for meeting
2) get a simple MPC example working?
3) email group to make sure they bring their laptops

# after
1) clean up website

==================
# High Priority #
==================
2) Make test functions
3) fix the solver options
4) fix integration of states for :ps methods

===================
# Medium Priority #
===================
17) make a better way to define diff eqs
18) add ! to allPlots()
19) add params to n.params
20) when integrating using :tm :trapezoidal method and minimizing the square of something it may ossilate
21) consider making dummy variables for state and control variables x1 x2, u1, u2 etc.
22) provide an initial guess
23) using the macro make two functions, one for JuMP and one for DiffEq
24) error checking in the DiffEq function
25) parse_DiffEq works OK except when in module
26) adding the NLexpression still does not work

=================
# Low Priority #
=================
1) warn user if final time constraint is active
6) try other linear solvers with Ipopt
12) allow user to select from using the IMatrix or quadrature
14) make functionality to easily compare all integration schemes
15) implement method with variable time steps for :tm methods, that is how it was before!
try to register functions with JuMP
17) change warning from Pkg.add(`` ...) to ""
18) allow for an variable array of dts
19) fix pgfplots in mpcdocs
20) fix OCPdef
