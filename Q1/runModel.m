
% Run the timestepping model
% Newton solve on each time step
% Input:
%   m - number of grid points on [0,1] space domain
%   tf - final time on the [0,tf] time domain
%   nt - number of time steps on the [0,tf] time domain

function runModel( m, tf, nt )

global epsilon alpha uL uR u0
global n h dt

% time and space discretisation
dt = tf/nt;
h = 1/(m-1);
n = m-2;
x = linspace(0,1,m)';

% boundary conditions
uL = 0;
uR = 1;

% PDE parameters
epsilon = 0.25;
alpha = 1;

% initial conditions
v0 = linspace(uL,uR,m)';
u0 = v0(2:m-1);
fprintf('Grid Size %d\n', n);
fprintf('Time %f\n',0);

plot(x,v0,'-ob'); hold on;

% time stepping
tic
for k = 1:nt
    t = k*dt;
    fprintf('Time %f\n',t);
    
%     [u{k},f,iter(k)] = newtonAlgorithm(@residual, u0, 1e-8, ...
%                         @fdJacobian, @linearSolve, 10);
    [u{k},f,iter(k)] = newtonAlgorithm(@residual, u0, 1e-8, ...
                        @tridiagonalJacobian, @sparseThomas, 10);
    u0 = u{k};
end
toc

disp(['Total number of Newton iterations: ' num2str(sum(iter))])

plot(x,[uL; u0; uR],'-r','LineWidth',2); 
legend('Initial solution u(x,0)', 'End solution u(x,t_f)','Location','SouthEast')
end


