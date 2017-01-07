Multiple States -> ex#3
==================================
In this example, we are going to approximate the `5th` order Taylor series polynomial for `sin()` and `cos()` expanded about `x=0`.

where:
 .. math:: sin(x) ≈ P_5(x) = x - \frac{x^3}{3!} + \frac{x^5}{5!}
 .. math:: cos(x) ≈ P_5(x) = 1 - \frac{x^2}{2!} + \frac{x^4}{4!}

and:
::

  t0 = Float64(0); tf = Float64(2);

Problem A
----------
with:
::

 ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=5,Nck=[200,3,100,5,100]);


.. image:: test3a.png


.. sidebar:: Comments on failing to calculate integral in Problem A

  This is expected. Looking at the smallest mesh grid size `Nck=1`, we can only expect to calculate a the integral for a `2*1-2=0th` order polynomial exactly.

Problem B
----------
with:
::

  ps, nlp = initialize_NLP(numStates=numStates,numControls=1,Ni=2,Nck=[4,4]);

.. image:: test3b.png
