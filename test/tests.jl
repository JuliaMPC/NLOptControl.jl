
function test1()
  x₀=-1;
  xₙ= 1;
  x = linspace(x₀,xₙ,100);
  α₁ =  3;
  α₂ = -3;
  α₃ = -8;
  α₄ =  7;
  γ = Poly([α₁,α₂,α₃,α₄]);
  y = polyval(γ,x)
  # construct polynomial approximation
  N = 2;
  τ = LGL_nodes(N);
  ω = LGL_weights(τ);
  D = LGL_Dmatrix(τ);
  # approximate the integral
  ∫γ = polyint(γ);
  Y = polyval(∫γ,x[end]) - polyval(∫γ,x[1]);
  C = Y - polyval(∫γ,x[end]); # find the constant of integration, using KNOWN solution from definite integral
  ∫y = polyval(∫γ,x) + C;
  fτ = polyval(γ,τ)
  ζ = cumsum(ω.*fτ,1)

  actual = ∫y[end]
  approx = ζ[end]
  return actual, approx
end
