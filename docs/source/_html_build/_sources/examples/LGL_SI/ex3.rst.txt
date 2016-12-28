Example 3
=========
In the third example, we approximate the derivative of a function by :ref:`diff_matrix`.


Test 3a
-------
In this test:

.. math:: y(x) = 3x^2-3x+3

.. image:: test3a.png

* We are able to determine the derivative exactly when they are linear functions is :math:`N = 3`

Test 3b
-------
In this test we increase the order of ``y(x)`` to:

.. math:: y(x) = 3x^3 + 3x^2-3x+3

.. image:: test3b.png

* When the derivative function becomes nonlinear

  * We can no longer calculate it exactly everywhere
  * There are only :math:`N_{t+1} = 4` node points
  * To calculate the derivative exactly we would need an :math:`∞` amount of :math:`N_{t+1}`

* We are still calculating the integral exactly and should be able to with :math:`N = 3` until :math:`x^5`

Test 3c
-------
In this test we increase the order of ``y(x)`` to:

.. math:: y(x) = 3x^4 + 3x^3 + 3x^2-3x+3

.. image:: test3c.png

* We are still calculating the integral exactly and should be able to with :math:`N = 3` until :math:`x^5`

Test 3d
-------
In this test we increase the order of ``y(x)`` to:

.. math:: y(x) =3x^5 + 3x^4 + 3x^3 + 3x^2-3x+3

.. image:: test3d.png

* We are still calculating the integral exactly with :math:`N = 3`!!
* The percent error is = 0.000000000000000000 %

Test 3e
-------
In this test we increase the order of ``y(x)`` to:

.. math:: y(x) =3x^6 + 3x^5 + 3x^4 + 3x^3 + 3x^2-3x+3

.. image:: test3e.png

* As expected, we are not still calculating the integral exactly with :math:`N = 3`!!
* The percent error is = -1.818181818181822340 %
