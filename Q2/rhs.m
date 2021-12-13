function g = rhs( x,y,z )
% Right hand side of our PDE
%   function g = rhs( x,y,z )
% Input: coordinate point (x,y,z)
% Output: evaluated function g

g = 2.0*exp( -0.5*( (x+3)^2 + (y-1)^2 + (z+4)^2 ) ) - 3.0*exp( -0.5*( (x-2)^2 + (y-1)^2 + (z-3)^2  ) );

end

