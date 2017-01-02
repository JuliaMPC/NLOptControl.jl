Neglecting Non-Collocated Point :math:`Y^{(k)}(τ)` -> ex#1
===========================================================
In this first example, we demonstrate the functionality using LGR nodes

where:
 .. math:: y(x) = -3x^5+3x^4+7x^3-8x^2-3x+3

Test 1a
-------
with:
  ::

    Nc = Int64(10); # number of collocation points in each interval
    Ni = Int64(4);  # number of intervals

.. image:: test1a.png


Test 1b
-------
with:
  ::

    Nc = Int64(3); # number of collocation points in each interval
    Ni = Int64(2);  # number of intervals

.. image:: test1b.png

Test 1c
-------
with:
  ::

    Nc = Int64(4); # number of collocation points in each interval
    Ni = Int64(2);  # number of intervals

.. image:: test1c.png

* Conclusions

  * It seems, that using the multiple interval formulation sacrifices the property where

      * We can calculate a ``Pth`` order polynomial exactly with :math:`2*N-2` collocation points
      * We do not approximate a ``5th`` order polynomial with ``6`` total collocation points
      * Looking at Test 1c, we can see that when we use ``N=4`` we calculate the intergral exactly

        * So the property applies only to each interval

  * Test 1b and Test 1c both show that we are not calculating the end point
