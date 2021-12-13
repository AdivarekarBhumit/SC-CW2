function [x,k] = jacobi(A, b, x0, tol, maxIter)

r = b - A*x0;
k=1;

% Split the matrix A = D + E to the diagonal and the rest
D = diag(diag(A));
E = A - D;

Bfun = @(X)(-D\(E*X));
N = size(A,1);
  
z = D \ b;
fprintf('Jacobi\nIter  RelResidual\n')
fprintf('------------------\n')

while (norm(r)/norm(x0) > tol && k < maxIter)
    x = Bfun(x0) + z;
    x0 = x;
    r = b - A*x0;
    fprintf('%3.0f   %1.3e\n',[k norm(r)/norm(x)])
    k=k+1;
end
