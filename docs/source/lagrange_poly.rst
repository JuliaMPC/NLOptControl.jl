.. _lagrange_poly:

Lagrange Interpolating Polynomials
==================================

Lagrange Polynomial Basics:

  * given :math:`(N+1)` unique data points

    * :math:`(x_0,y_0),(x_1,y_1),....,(x_N,y_N)`
  * we can create an :math:`N^{th}` order Lagrange interpolating polynomial

.. math:: P_n(x) = \sum_{i=0}^N \mathcal{L}_i(x)f(x_i)

where,
  .. math::
      :nowrap:

      \begin{eqnarray}
        f(x_0) = y_0\\
        f(x_1) = y_1\\
        .\\
        .\\
        f(x_i) = y_i\\
        .\\
        f(x_N) = y_N
       \end{eqnarray}

So, we are just multiplying by the given :math:`y_i` values.

Lagrange Polynomial Basis Polynomials
---------------------------------------
.. math:: \mathcal{L}_i(x)=\prod_{\substack{j=0 \\ j\neq i}}^{N}\frac{x-x_j}{x_i-x_j}

so expanding this,
  .. math::
      :nowrap:

      \begin{eqnarray}
      \mathcal{L}_i(x) &=\frac{x-x_0}{x_i-x_0}\frac{x-x_1}{x_i-x_1}...\\
                       &...\frac{x-x_{i-1}}{x_i-x_{i-1}}...\\
                       &...\frac{x-x_{i+1}}{x_i-x_{i+1}}...\\
                       &...\frac{x-x_N}{x_i-x_N}
      \end{eqnarray}

Notice that we do not include the term where :math:`i==j`!


These equations where turned into a function:
::

  """
  L = lagrange_basis_poly{T<:Number}(x::Float64,x_data::AbstractArray{T},idx::Int64,N::Int64)
  --------------------------------------------------------------------\n
  Author: Huckleberry Febbo, Graduate Student, University of Michigan
  Date Create: 12/26/2016, Last Modified: 12/26/2016
  --------------------------------------------------------------------\n
  # Input Arguments
  * `x::Float64`: point to approximate function value at
  * `x_data::Array{Float64}`: x data to used calculate basis polynomials
  * `N::Int64`: order of Lagrange interpolating polynomial

  # Output Arguments
  * `L::Array{Float64}`: Lagrange basis polynomials

  A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

  """

The development of this function can be seen here:
::
  notebook(dir=Pkg.dir("NLOptControl/examples/LIP/lagrange_basis_ploy_dev"))


More information `here <https://en.wikipedia.org/wiki/Vandermonde_matrix>`_
