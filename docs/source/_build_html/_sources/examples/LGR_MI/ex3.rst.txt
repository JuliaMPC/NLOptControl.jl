Approximation of State Derivative at of Mesh Grids -> ex#3
==========================================================
In this example, the state derivative at the end of each mesh.

where:
 .. math:: y(x) = -0.3x^2+3x-0.3

with:
  ::

    Nc = Int64(3); # number of collocation points in each interval
    Ni = Int64(2);  # number of intervals


.. image:: test3a.png


.. sidebar::  Why Do We Need This Derivative At the End Point?

 Actually, we do not need this. There is no constraint on state dynamics at the end of each mesh grid using the method discussed in :cite:`d-garg2011advances`.  This is an important point, that is described here :ref:`hp_description`.


Now, we will look at the :math:`\mathbf{D}` matrix used to calculate the derivatives above:
::

     D =
     4×4×2 Array{Float64,3}:
     [:, :, 1] =
      -1.0         1.50639   -1.10639    0.6
      -0.210639   -0.155051   0.713568  -0.347878
       0.0506395  -0.233568  -0.644949   0.827878
      -0.0666667   0.276429  -2.00976    1.8

     [:, :, 2] =
      -1.0         1.50639   -1.10639    0.6
      -0.210639   -0.155051   0.713568  -0.347878
       0.0506395  -0.233568  -0.644949   0.827878
      -0.0666667   0.276429  -2.00976    1.8

Notice that for each interval the :math:`\mathbf{D}` matrix is actually identical. This is quite an interesting observation indeed, because different inputs where used to calculate it, these are the nodes.

For the first interval:
::

   0.0
   1.77526
   4.22474
   5.0

For the second interval:
::

  5.0
  6.77526
  9.22474
  10.0


These nodes depend on the interval :math:`t_0->t_f` as well as the :math:`\tau`:
::

 -1.0
 -0.289898
  0.689898

Which are the LGR nodes when :math:`N_c=3`

So, it seems that maybe we can calculate the weights beforehand as well as the :math:`\mathbf{D}` matrix and cache the result.

Neglecting Derivative At End Of Mesh
---------------------------------------
For the purposes of using this method for control, again we do not need to calculate the derivative of the state at the ends of each mesh. So, we can remove the bottom row of the :math:`\mathbf{D}` matrix as:
::

  D =
  [:, :, 1] =
   -1.0         1.50639   -1.10639    0.6
   -0.210639   -0.155051   0.713568  -0.347878
    0.0506395  -0.233568  -0.644949   0.827878

  [:, :, 2] =
   -1.0         1.50639   -1.10639    0.6
   -0.210639   -0.155051   0.713568  -0.347878
    0.0506395  -0.233568  -0.644949   0.827878

.. image:: test3b.png

So, at the end of each mesh grid, we still approximate the state, but neglect it's derivative.


.. rubric:: References

.. bibliography:: zref.bib
   :labelprefix: D
   :keyprefix: d-
