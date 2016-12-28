Example 1
=========
In the first example, we borrow a problem from `Wikipedia <https://en.wikipedia.org/wiki/Gaussian_quadrature>`_.

where:
 .. math:: y(x) = 7x^3-8x^2-3x+3


.. image:: test1.png


.. sidebar:: Difference between the Wikipedia Example and this Example

  The difference between Wikipedia example and this one is that the Wikipedia example uses Gauss-Legendre Quadrature while the code developed in this package uses Legendre-Pseudospectral Method with Lagrange-Gauss-Lobatto (LGL) nodes. Information on the difference between these methods and many more can be found in both :cite:`b-Shen2011` and :cite:`b-herber2015basic`.


* Conclusions

  * We are able to exactly determine the integral

.. rubric:: References

.. bibliography:: zref.bib
    :labelprefix: B
    :keyprefix: b-
