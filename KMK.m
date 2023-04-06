function y = KMK(t,x,beta,gamma)

% INPUT: t is a a real value indicating time
% x is a column vector of size 3 x 1
% beta, gamma are parameters of the KMK equations

% OUPUT: y is a column vector of size 3 x 1 that gives
% the right hand side of the KMK equations

S = x(1,1);
I = x(1,2);
R = x(1,3);

y = [-beta*S*I, beta*S*I - gamma*I, gamma*I];

end
