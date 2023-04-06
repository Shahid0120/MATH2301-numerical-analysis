function [tout,yout] = EulerSolver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE dy/dt = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a; b; c] is a column vector that contain the
% initial values y(0)= y0 = [a;b;c].
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

tout = t; % ensures tout is a column vector 
y0 = y0(:); % ensures y0 is a column vector
yout = zeros(length(y0), length(t)); % initialize the solution matrix
yout(:, 1) = y0; % set the initial values

for i = 1:length(t)-1
    h = t(i+1) - t(i);
    yout(:, i+1) = yout(:, i) + h*f(t(i), yout(:, i)); % Euler's method
end

end
