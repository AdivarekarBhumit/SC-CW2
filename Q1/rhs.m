
% Right hand side of the PDE system udot = f(u)
% PDE: du/dt = epsilon d2u/dx2 + alpha (du/dx)^2

function r = rhs( u )

global n h epsilon alpha uL uR

v = [ uL; u; uR ];

r = zeros(n+2,1);

for i = 2:n+1
    r(i) = ( epsilon*(v(i+1)-2*v(i)+v(i-1)) + 0.25*alpha*(v(i+1)-v(i-1))^2 )/h^2;
end

r = r(2:n+1);

end
