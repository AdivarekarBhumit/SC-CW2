function [ x,y,z, u ] = runFDM( mIn )
% Solve our nonlinear 3D FDM model 
%   function [ x,y, u ] = runFDM( mIn )
% Input: mIn - grid size in each coordinate direction
% Output: x,y,z,u - 3d grid (x,y) and solution data

global m N h x y z
global A loc nf

close('all');

m = mIn;
[x,y,z] = meshgrid( linspace(-10,10,m ) ); % 3d mesh for entire domain
h = 20/(m-1);
N = (m-2)^3; % unknown d.o.f

fprintf('\n Grid size %d x %d x %d\n',m,m,m);
fprintf(' Nonlinear system size %d\n\n', N );

[ A, nz ] = sparsePattern3DFD( m ); % 3d FDM sparsity pattern %
%Uncommented for Q2 b
[ loc, nf ] = sparseJac( A,N );     % sparse function calls for J evaluation
%Uncommented for Q2 b

u0 = zeros( N,1 ); % initial guess for interior solution

% call Newton solver: 
%  u0 - initial guess
%  @fdm3D - the function defining the nonlinear equations
%  @fdJacobian - a function that computes the Jacobian matrix
%  @linearSolve - a function that solves the linear system J.d=-f
%[ u,f,iter ] = newtonAlgorithm( @fdm3D, u0, 1e-8, @buildJacobian, @linearSolve, 40 );
[ u,f,iter, J ] = newtonAlgorithmJ( @fdm3D, u0, 1e-8, @buildJacobian, @linearSolve, 40 );

b = ones(N, 1);
tic;
for i=1:100
   directLU(J,b);
end
tA1o = toc;

fprintf('Direct LU Time %f\n', tA1o);

tic;
for i=1:100
   x1gs = jacobi(J, b, b, 1e-6, 1000);
end
tA1o = toc;
fprintf('Jacobi Time %f\n', tA1o);

tic;
for i=1:100
   x1gs = gaussseidel(J, b, b, 1e-6, 1000);
end
tA1o = toc;
fprintf('GaussSidel Time %f\n', tA1o);

end
