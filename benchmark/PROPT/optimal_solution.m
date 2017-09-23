function [h,v,u] = optimal_solution(t,v0,h0,ts)
% actual optimal solution
  if t < ts
    h = -3/4*t^2 + v0*t + h0;
    v = -3/2*t + v0;
    u = 0;
  else
    h = 3/4*t^2 + (-3*ts + v0)*t + 3/2*ts^2 + h0;
    v = 3/2*t + (-3*ts + v0);
    u = 3;
  end
end

