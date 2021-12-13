
% Nonlinear residual for the PDE: udot = rhs(u)
% PDE: du/dt = epsilon d2u/dx2 + alpha (du/dx)^2
% Implicit Euler in time 

function f = residual( u )

global dt u0 

r = rhs( u );

f = (u-u0)/dt - r;

end

