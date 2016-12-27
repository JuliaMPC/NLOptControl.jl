.. _lagrange_poly_fun:

Functionality
*************
The basic description of this functionality is detailed here :ref:`lagrange_poly`

`lagrange_basis_poly()`
=======================

The Lagrange basis polynomial equations where turned into a function.


`interpolate_lagrange()`
========================

The interpolation functionality was pushed to a lower level. This allows the user to easily use code to interpolate a polynomial.

The **development** of these function can be:

  *  Viewed remotely on  using the `jupyter nbviewer <http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/LIP/lagrange_basis_poly_dev.ipynb>`_.
  *  Viewed locally and interacted using IJulia

    To do this in julia type:
    ::

      using IJulia
      notebook(dir=Pkg.dir("NLOptControl/examples/LIP/lagrange_basis_poly_dev"))

Examples
*********
.. toctree::
   :maxdepth: 1

   ex1
   ex2
   ex3

These examples can be:

  *  Viewed remotely on  using `the jupyter nbviewer <http://nbviewer.jupyter.org/github/huckl3b3rry87/NLOptControl.jl/blob/master/examples/LIP>`_.
  *  Viewed locally and interacted using IJulia

  To do this in julia type:
  ::

    using IJulia
    notebook(dir=Pkg.dir("NLOptControl/examples/LIP/"))
