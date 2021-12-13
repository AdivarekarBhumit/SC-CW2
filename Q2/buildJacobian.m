
function J =  buildJacobian( n, u0, f0, fnon )

% Numerical Jacobian for sparse matrices
% Sparsity pattern precomputed and stored in (global) A
% Column groupings precomputed by sparseJac.m in (global) loc,nf
%
% function J =  fd2dJacobian( u )
%  Input: n - number of equations
%         u0 - base solution values
%         f0 - base function values
%         fnon - nonlinear residual function
%  Output: J - Jacobian matrix in sparse form

global loc nf A

J = A;

he = 10*sqrt(eps);

count = 1;

for kf = 1:nf   % function evaluation

    u = u0; % base state

    Jstart = count;
    while( count<n+1 && loc(count,1)==kf )   % next variable
        j = loc(count,2);
        u(j) = u(j) + he; % perturbation
        count = count + 1;
    end
    Jend = count - 1;

    f = feval( fnon, u ); % function call
    jac = ( f - f0 )/he;   % difference approximation

    for k = Jstart:Jend   % variable
        j = loc(k,2);
        [iloc,jloc] = find(A(:,j));
        for m = 1:length(iloc)
            J(iloc(m),j) = jac(iloc(m));
        end
    end
    
end

end
