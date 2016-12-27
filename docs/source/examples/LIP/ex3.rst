Runge's Phenomena
=================

This example investigates `Runge's phenomena <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.

where:
 .. math:: y(x) = \frac{1}{1+25x^2}

and the interval from ``x0=-1`` to ``xf=1``

with:
  ::

    N = order of Lagrange Polynomial
    x_data = linspace(x0,xf,N+1);

.. image:: test3a.png


* Conclusions

  * Be careful not to use too high of an ``N``

.. sidebar:: To Mitigate Runge's Phenomenon

  * Could sample more near the end points
  * Could use `Chebyshev nodes <https://en.wikipedia.org/wiki/Chebyshev_nodes>`_
