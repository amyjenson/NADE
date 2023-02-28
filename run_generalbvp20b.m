close all
clear all

% This script run the function generalbvp.m which 
% approximates the solution to the boundary value problem 
% u''(x) + p(x) u'(x) + q(x) u(x) = f(x) , u(xL) = alpha, u(xR) = beta
% using a centered finite difference method (see Problem 11 in A2).
%
%function [x,U] = bvpq(m, xL,xR, p, q, f, alpha, beta) 


f = @(x) 0*x;
xL = 0;
xR = 1;
alpha = 1;
beta = 0;
p =@(x) -20;
q =@(x) 0*x; 

for i = [2500] %[3,5,10,20,50,200,1000]
     [x,U] = generalbvp(i, xL, xR, p, q, f, alpha, beta);
     plot(x,U)
     hold on
end   

u_ex =@(x) 1-((1-exp(20*x))/(1-exp(20)));
u_ex = u_ex(x);

plot(x, u_ex)

legend('m=3', 'm=5', 'm=10', 'm=20', 'm=50', 'm=200','m=1000', 'exact solution' , 'Location', 'Southwest')

xlabel('x')
ylabel('u')
title('Numerical solutions for different values of m plotted with exact solution')