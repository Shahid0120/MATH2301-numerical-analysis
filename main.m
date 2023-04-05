beta = 0.6;
gamma = 1/3;

f = @(t,x) KMK(t,x,beta,gamma);

%Parameters can be changed according to Task 6
k = 2; 
h = 10^(-k);
tfinal = 150;
t = 0:h:tfinal;
y0 = [1;1.27*1e-6;0];

solve = input('What method do you want to use? Euler, RK4 or IRK4. ', 's');
if strncmpi(solve,'Euler',5)
    [tout, Y] = EulerSolver(f,t,y0);
elseif strncmpi(solve, 'RK4',3)
    [tout , Y] = RK4Solver(f,t,y0);
elseif strncmpi(solve,'IRK4',4)
    [tout , Y] = IRK4Solver(f,t,y0);
else
    fprintf('You did not select a method available, try again.\n')
    return
end

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);
[tmout , Ym] = ode45(f,t,y0,options);

figure(1)
subplot(2,1,1)
plot(tmout,Ym)
xlabel('t');
legend('S', 'I', 'R');
title('ode45');
subplot(2,1,2)
plot(tout,Y)
xlabel('t');
legend('S', 'I', 'R');
title(solve);


err = max(max(abs(Y-Ym')));
fprintf('\nh = %1.1e \t error = %.10e\n\n', h, err);
