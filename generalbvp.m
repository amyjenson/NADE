function [x,U] = generalbvp(m, xL,xR, p, q, f, alpha, beta) 
% This script approximates the solution to the boundary value problem 
% u''(x) + p u'(x) + q u(x) = f(x) , u(xL) = alpha, u(xR) = beta
% using a centered finite difference method (see Problem 11 in A2).
%
%function [x,U] = generalbvp(m, xL,xR, p, q, f, alpha, beta) 
%
% Example:  Attempt to solve m=20 test case problem.
%   >> f = @(x) 2 - 12 * x + 12 * x.^2;
%   >> [x,U] = generalbvp(20,0, 1, p, q, f , 1,0 );

% set up grid
h = (xR-xL)/(m+1);     % length of step
x = xL:h:xR;     % length of interval is m+2 

% assemble linear system  A U = F
A =  -2/h^2 * eye(m);
for i=1:m
    A(i,i) = A(i,i) + q(x(i));
end
for j = 1:m-1
    A(j,j+1) = 1/h^2 - p(x(j))/2*h;
    A(j+1,j) = 1/h^2 + p(x(j+1))/2*h;
end
% evaluate f at *all* grid points, then correct values
F = f(x(2:m+1))'; 
F(1) = F(1) - (1/h^2 - p(x(1))/2*h)*alpha;
F(m) = F(m) - (1/h^2 + p(x(m))/2*h)*beta;
% solve the linear system
U = A\F;           % numerical solution at interior points
U = [alpha U' beta]; % whole solution including boundary vals
