**********************************
NLP Problem Initialization
**********************************

Here we are developing software to set up and keep track of all of the variables in the problem.

Functionality
#############

`ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=1,Nck=[2])`
*******************************************************************
* To initialize the problem

To use:
::

  using NLOptControl
  ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=2,Nck=[2,3]);

ps = pseudospectral method related data:
::

  NLOptControl.PS_data
    Nck: [2,3]
    Ni: 2
    τ: Array{Float64,1}[[-1.0,0.333333],[-1.0,-0.289898,0.689898]]
    ts: Array{Float64,1}[[-0.0,0.0],[-0.0,-0.0,0.0]]
    ω: Array{Float64,1}[[0.5,1.5],[0.222222,1.02497,0.752806]]
    ωₛ: Array{Float64,1}[[0.0,0.0],[0.0,0.0,0.0]]
    t0: 0.0
    tf: 0.0
    DMatrix: Array{Float64,2}[[0.0 0.0 0.0; 0.0 0.0 0.0],[0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]]
    IMatrix: Array{Float64,2}[[0.0 0.0; 0.0 0.0],[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]]
    stateMatrix: Array{Float64,2}[[0.0; 0.0; 0.0],[0.0; 0.0; 0.0; 0.0]]
    controlMatrix: Array{Float64,2}[[0.0; 0.0],[0.0; 0.0; 0.0]]

nlp = nonlinear programming problem related data:
::

  NLOptControl.NLP_data
  numStates: 1
  numControls: 1
  numStatePoints: [3,4]
  numControlPoints: [2,3]
  lengthStateVector: 7
  lengthControlVector: 5
  lengthDecVector: 14
  timeStartIdx: 13
  timeStopIdx: 14
  stateIdx: Tuple{Int64,Int64}[(1,3),(4,7)]
  controlIdx: Tuple{Int64,Int64}[(8,9),(10,12)]
  stateIdx_all: Tuple{Int64,Int64}[(-99,-99)]
  controlIdx_all: Tuple{Int64,Int64}[(-99,-99)]
  stateIdx_st: Tuple{Int64,Int64}[(-99,-99)]
  controlIdx_ctr: Tuple{Int64,Int64}[(-99,-99)]
  decisionVector: [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

Comments
---------

  * Lots of zeros right now because we still have a lot of initializing to do.
  * `numControls` and `numStates` can both be greater than `1`.
  * The value of `Ni` must match the length of `Nck`


`generate_Fake_data(nlp,ps,γ)`
*******************************

* Generating some data to play with is useful:

`nlp2ocp(nlp,ps)`
******************
* Turning the nonlinear-programming problem into an optimal control problem

  * This function basically takes a large design variable an sorts it back into the original variables

`ζ, approx_int_st = integrate_state(ps,nlp)`
*********************************************
Approximates integral.

To use this function  with Gaussian Quadrature (the default method):
::

  ζ, approx_int_st = integrate_state(ps,nlp)


To use this function  with the LGRIM:
::

  integrate_state(ps,nlp;(:mode=>:LGRIM))


`dζ = differentiate_state(ps,nlp)`
***********************************
Approximate derivatives.

* Currently only using LGRDM as a method.


Examples
#########
In these examples we use:
  * Demonstrate functionality to setup optimal control problem
  * Also include the development scripts of these functions
  * There is not a webpage for all examples, but the interested user can check out them out using IJulia

.. toctree::
   :maxdepth: 2

   ex2

These examples can be:
  *  Viewed remotely on  using the `jupyter nbviewer <http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/NLP2OCP>`_.
  *  Viewed locally and interacted using IJulia

      To do this in julia type:
      ::

        using IJulia
        notebook(dir=Pkg.dir("NLOptControl/examples/NLP2OCP/"))


.. sidebar::  Researchers at the University of Florida

  Describe this method in many papers including :cite:`d-darby2011hp,d-garg2011advances,d-garg2010unified,d-garg2009overview`.


.. rubric:: References

.. bibliography:: zref.bib
    :labelprefix: D
    :keyprefix: d-
