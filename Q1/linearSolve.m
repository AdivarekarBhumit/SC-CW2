function x = linearSolve( A, b )

% function x = linearSolve( A, b )
%  Linear solution algorithm for the
%  system Ax=b - dense linear algebra
% Input: nxn dense matrix A
%        n-dim column vector b
% Output: n-dim column vector x

n = length(b);

for i = 1:n-1
  r = 1/A(i,i);
  for j = i+1:n
    d = r*A(j,i);
    for k = i:n
      A(j,k) = A(j,k) - d*A(i,k);
    end
    b(j) = b(j) - d*b(i);
  end
end 

x = zeros(n,1);

x(n) = b(n)/A(n,n);
for i = n-1:-1:1
  x(i) = ( b(i) - A(i,i+1:n)*x(i+1:n) ) / A(i,i);
end   

end