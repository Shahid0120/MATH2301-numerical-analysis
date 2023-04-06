function [tout,yout] = IRK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a; b; c] is a column vector that contain the
% initial values y(0)=y0=[a;b;c].
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

tout = t; % ensures tout is a column vector 
y0 = y0(:); % ensures y0 is a column vector
yout = zeros(length(y0), length(t)); % initialize the solution matrix
yout(:, 1) = y0; % set the initial values

% Need to define the IRK Tableau (Two-Stage of Order 4)

c = [1/2 - sqrt(3)/6, 1/2 + sqrt(3)/6];
b = [1/2, 1/2];
A = [1/4, 1/4 - sqrt(3)/6; 1/4 + sqrt(3)/6, 1/4];


%For Loop for determining the iteration solution

for k = 1:(length(t) - 1)
    
    h = t(k+1)-t(k);

    XI1 = @(xi1, xi2) yout(:, k) + h*A(1, 1)*f(t(k) + c(1)*h, xi1) + h*A(1, 2)*f(t(k) + c(2)*h, xi2) - xi1;
    XI2 = @(xi1, xi2) yout(:, k) + h*A(2, 1)*f(t(k) + c(1)*h, xi1) + h*A(2, 2)*f(t(k) + c(2)*h, xi2) - xi2;
    
    p = @(p) [XI1(p(:,1),p(:,2)),XI2(p(:,1),p(:,2))];
    
    x0 = [yout(:, k), (yout(:, k) + h*f(t(k), yout(:, k)))];
    
    options=optimset('Display','off','TolFun',1e-32,'TolX',1e-32);
    
    psi = fsolve(p,x0,options); 
   
    yout(:, k + 1) = yout(:, k) + h*b(1)*f(t(k) + c(1)*h, psi(:,1)) + h*b(2)*f(t(k) + c(2)*h, psi(:,2));
end

end
