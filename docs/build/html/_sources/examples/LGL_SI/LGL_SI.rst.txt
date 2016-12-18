Single Interval
===============

Basic Problem Definition
-------------------------
The code developed in this package uses the Legendre-Pseudospectral Method with Lagrange-Gauss-Lobatto (LGL) nodes. A basic description of this implementation presented in this documentation at :ref:`pseudospectral` and more much more detailed information can be found in both :cite:`Shen2011` and :cite:`herber2015basic`.

Examples
--------
In these examples we use:
  * Legendre-Gauss-Lobatto (LGL) nodes
  * Single interval approximations
  * Approximate integrals in the range of ``[-1,1]``
  * Approximate derivatives in the range of ``[-1,1]``

These examples can be:
  *  Viewed remotely on  using the `jupyter nbviewer <http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/LGL_SI>`_.
  *  Viewed locally and interacted using IJulia

      To do this:
      ::

         julia>using IJulia
         julia>notebook(dir=Pkg.dir("NLOptControl/examples/LGL_SI/"))


.. toctree::
   :maxdepth: 2

   ex1
   ex2
   ex3
