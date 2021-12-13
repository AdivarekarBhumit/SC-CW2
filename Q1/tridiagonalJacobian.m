
function J =  tridiagonalJacobian( n, u0, f0, fnon )

% Numerical tridiagonal Jacobian
%  function J =  tridiagonalJacobian( u )
% Input: n - number of equations
%        u0 - base solution values
%        f0 - base function values
%        fnon - nonlinear residual function
% Output: J - Jacobian matrix in sparse tridiagonal form

nz = 3*n-2;
irow = zeros(nz,1);
jcol = zeros(nz,1);
a = zeros(nz,1);

del = 10*sqrt(eps);

z = 0;

for k = 1:3   % function evaluation

    u = u0; % base state

    for j = k:3:n   % variable
        u(j) = u(j) + del; % perturbation
    end

    f = feval( fnon, u ); % function call
    jac = ( f - f0 )/del;   % difference approximation

    for j = k:3:n   % variable
        for i = max(1,j-1):min(n,j+1)   % equation
            z = z+1;
            irow(z) = i;
            jcol(z) = j;
            a(z) = jac(i);
        end
    end    
end

J = sparse(irow,jcol,a,n,n,nz);

end
