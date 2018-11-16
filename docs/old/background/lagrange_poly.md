Lagrange Interpolating Polynomials
==================================

Definition
----------

> -   given $(N+1)$ unique data points
>
>     \* $(x_0,y_0),(x_1,y_1),....,(x_N,y_N)$
> -   we can create an $N^{th}$ order Lagrange interpolating polynomial

$$P_n(x) = \sum_{i=0}^N \mathcal{L}_i(x)f(x_i)$$

where,
:   $$:nowrap:$$$$\begin{eqnarray}
      f(x_0) = y_0\\
      f(x_1) = y_1\\
      .\\
      .\\
      f(x_i) = y_i\\
      .\\
      f(x_N) = y_N
     \end{eqnarray}$$

So, we are just multiplying by the given $y_i$ values.

Lagrange Basis Polynomials
--------------------------

More information on Lagrange Basis Polynomials is
[here](https://en.wikipedia.org/wiki/Vandermonde_matrix)

$$\mathcal{L}_i(x)=\prod_{\substack{j=0 \\ j\neq i}}^{N}\frac{x-x_j}{x_i-x_j}$$

so expanding this,
:   $$:nowrap:$$$$\begin{eqnarray}
    \mathcal{L}_i(x) &=\frac{x-x_0}{x_i-x_0}\frac{x-x_1}{x_i-x_1}...\\
                     &...\frac{x-x_{i-1}}{x_i-x_{i-1}}...\\
                     &...\frac{x-x_{i+1}}{x_i-x_{i+1}}...\\
                     &...\frac{x-x_N}{x_i-x_N}
    \end{eqnarray}$$

Notice that we do not include the term where $i==j$!

Please see lpf for details on implementation.
