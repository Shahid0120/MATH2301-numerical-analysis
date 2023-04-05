function [tout,yout] = RK4Solver(f,t,y0)
% INPUT: f(t,y) is an anonymous function that defines
%        the right-hand side of the ODE ydot = f(t,y)
%        t =[t0 t1 ... tfinal] is a vector of grid points
%        with length N
%        y0=[a; b; c] is a column vector that contain the
%        initial values y(0)=y0 = [a;b;c].
% OUTPUT:tout is a column vector of grid points.
%        yout is an 3 x N matrix containing the solution
%        at different grid points.

N = length(t);
y0 = y0(:);
yout = zeros(length(y0),N);
yout(1,:) = y0;
tout = t(:);
h = t(2) - t(1);

for i=1:N-1
    k1 = f( t(i), yout(:,i));
    k2 = f( t(i) + 0.5*h, yout(:,i) + 0.5*h*k1 );
    k3 = f( t(i) + 0.5*h, yout(:,i) + 0.5*h*k2 );
    k4 = f( t(i) + h, yout(:,i) + h*k3 );
    yout(:,i+1) = yout(:,1) + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end

end 
