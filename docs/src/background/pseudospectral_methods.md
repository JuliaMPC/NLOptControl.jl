Pseudospectral Methods
======================

Change of Interval
------------------

To can change the limits of the integration (in order to apply
Quadrature), we introduce $\tau \in [-1,+1]$ as a new independent
variable and perform a change of variable for $t$ in terms of $\tau$, by
defining:

> $$\tau = \frac{2}{t_{{N}_{t}}-t_0}t - \frac{t_{N_t}+t_0}{t_{N_t}-t_0}$$

Polynomial Interpolation
------------------------

Select a set of $N_t+1$ node points:

> $$\mathbf{\tau} = [\tau_0,\tau_1,\tau_2,.....,\tau_{N_t}]$$

-   These none points are just numbers
    -   Increasing and distinct numbers $\in [-1,+1]$

A *unique* polynomial $P(\tau)$ exists (i.e. $\exists! P(\tau)$) of a
maximum degree of $N_t$ where:

> $$f(\tau_k)=P(\tau_k),\;\;\;k={0,1,2,....N_t}$$

-   So, the function evaluated at $\tau_k$ is equivalent the this
    polynomial evaluated at that point.

But, between the intervals, we must approximate $f(\tau)$ as:

> $$f(\tau) \approx P(\tau)= \sum_{k=0}^{N_t}f(\tau_k)\phi_k(\tau)$$

with $\phi_k(•)$ are basis polynomials that are built by interpolating
$f(\tau)$ at the node points.

Approximating Derivatives
-------------------------

We can also approximate the derivative of a function $f(\tau)$ as:

$$\frac{\mathrm{d}f(\tau)}{\mathrm{d}\tau}=\dot{f}(\tau_k)\approx\dot{P}(\tau_k)=\sum_{i=0}^{N_t}D_{ki}f(\tau_i)$$

With $\mathbf{D}$ is a $(N_t+1)\times(N_t+1)$ differentiation matrix
that depends on:

> -   values of $\tau$
> -   type of interpolating polynomial

Now we have an approximation of $\dot{f}(\tau_k)$ that depends only on
$f(\tau)$!

Approximating Integrals
-----------------------

The integral we are interested in evaluating is:

$$\int_{t_0}^{t_{N_t}}f(t)\mathrm{d}t=\frac{t_{N_t}-t_0}{2}\int_{-1}^{1}f(\tau_k)\mathrm{d}\tau$$

This can be approximated using quadrature:

$$\int_{-1}^{1}f(\tau_k)\mathrm{d}\tau\sum_{k=0}^{N_t}w_kf(\tau_k)$$

where $w_k$ are quadrature weights and depend only on:

> -   values of $\tau$
> -   type of interpolating polynomial

Legendre Pseudospectral Method
------------------------------

-   Polynomial

Define an N order Legendre polynomial as:

> $$L_N(\tau) = \frac{1}{2^NN!}\frac{\mathrm{d}^n}{\mathrm{d}\tau^N}(\tau^2-1)^N$$

-   Nodes

$$:nowrap:$$$$\begin{equation}
  \tau_k = \left \{
  \begin{aligned}
    &-1, && \text{if}\ k=0 \\
    &\text{kth}\;\text{root}\;of \dot{L}_{N_t}(\tau), && \text{if}\ k = {1, 2, 3, .. N_t-1}\\
    &+1\, && \text{if}\ k = N_t
  \end{aligned} \right.
\end{equation}$$

-   Differentiation Matrix
-   Interpolating Polynomial Function

