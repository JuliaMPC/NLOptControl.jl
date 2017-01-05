NLP and OCP Functionality -> ex#2
==================================
In this example, we are continuing to work on preparing the code for use with optimization by creating higher level functionality. Examine the IJulia notebook to see code.

Aside from using the new data structures we demonstrate:

* Using the integration matrix for the first time
* Using the higher level functionality

where:
 .. math:: y(x) = -0.3x^2+3x-0.3

Single Interval Problem
------------------------
with:
::

  ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=1,Nck=[2]);

.. image:: test2a.png


Multiple Interval Problem A
----------------------------
with:
::

  ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=3,Nck=[2,3,5]);

.. image:: test2b.png


Multiple Interval Problem B
----------------------------
with:
::

  ps, nlp = initialize_NLP(numStates=1,numControls=1,Ni=5,Nck=[200,3,100,5,100]);


.. image:: test2c.png
