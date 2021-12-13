function [ x,f,normf ] = lineSearch( fnon, x0, normf0, delta )

% Line search algorithm for Newton update step
%  function [ x,f,normf ] = lineSearch( fnon, x0, normf0, delta )
% Input: fnon - function handle for nonlinear system
%        x0 - solution at start of Newton step
%        normf0 - norm of nonlinear system evaluated at x0
%        delta - update vector computed by Newton algorithm
% Output: x - final updated point
%         f - final function value
%         normf - final norm of nonlinear system

maxStep = 4; % max number of line search steps per iteration

step = 0;
lambda = 1;
x1 = x0 + delta;        % first update of x
f1 = feval(fnon,x1);   % new function values
normf1 = norm(f1);

fprintf('     alpha       |f(x)|\n')
fprintf('   %8.4f %12.6f\n',lambda,normf1)

while ( normf1 > normf0 )  && ( step < maxStep ) % line search
  step = step + 1;
  lambda = lambda/2;
  x1 = x0 + lambda*delta;  % update x
  f1 = feval(fnon,x1);
  normf1 = norm(f1);
  
  fprintf('   %8.4f %12.6f\n',lambda,normf1)
end

if( ( step == maxStep ) && ( normf1 > normf0 ) )
    fprintf('   Reached maxSteps on this Newton step\n');
end

x = x1;   % accept final value
f = f1;
normf = normf1;

end
