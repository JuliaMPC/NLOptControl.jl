**********************
LGR Multiple Interval
**********************

Now, using Legendre-Gauss-Radau (LGR) points with multiple intervals to calculate the integral and the derivative of a known polynomial function. This example demonstrates using the Legendre-Gauss-Radau (LGR) points to calculate the integral and the derivative of a known polynomial function using a **multiple interval approach**.


.. sidebar::  Researchers at the University of Florida

  describe this method in many papers including :cite:`c-darby2011hp,c-garg2011advances,c-garg2010unified,c-garg2009overview`.

Functionality
#############

`asa()`
********

Examples
#########
In these examples we use:
  * Legendre-Gauss-Radau (LGR) nodes
  * Multiple interval approximations
  * Approximate integrals in the range of ``[x0,xf]``
  * Approximate derivatives in the range of ``[x0,xf]``

.. toctree::
   :maxdepth: 2

   ex1
   ex2
   ex3
   ex4

These examples can be:
  *  Viewed remotely on  using the `jupyter nbviewer <http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/LGR_MI>`_.
  *  Viewed locally and interacted using IJulia

      To do this in julia type:
      ::

        using IJulia
        notebook(dir=Pkg.dir("NLOptControl/examples/LGR_MI/"))

.. rubric:: References

.. bibliography:: zref.bib
    :labelprefix: C
    :keyprefix: c-
