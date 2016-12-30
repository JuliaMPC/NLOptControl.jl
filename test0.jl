using ForwardDiff
using NLOptControl

# Derivative of Each Element
Nc = 2; x_data = [1,2,3]; j=1; x = 1;
d = ForwardDiff.derivative(lagrange_basis_poly, x, x_data, Nc, j)
