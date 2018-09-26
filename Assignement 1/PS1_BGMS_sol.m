%Solution of PSET 1 for
%Lionel Brodard, Tomas Giro, Antoine Marchal, Ivan Sch�nenberger

%Exercise 1
%Declaration of variables as stated in the exercise
S=100;
K=100;
T=1.5;
r=0.04;
%Convergence threshold
eps=0.01;
%First guess x_0 for Newton's algorithm
sigma=1;
%While convergence is not reached, keep iterating
while abs(f(S,K,T,sigma,r)-10.78)>eps
    sigma=iterateStep(S,K,T,sigma,r);
end
fprintf('Exercise 1 result:\n')
fprintf('The implied volatility is %.2f \n',sigma)
fprintf('---------------\n')

%Exercise 2

x_0=[3;3];
[x,n]=newton(x_0,1e-6);
fprintf('Exercise 2 part a) result:\n')
fprintf('The minimum of the function (x1-2)^4+(x1-2*x2)^2 is (%.4f;%.4f) \n',x(1),x(2))
fprintf('%d iterations were needed for Newtons method to converge with starting point (%.2f;%.2f)\n',n,x_0(1),x_0(1))
fprintf('---------------\n')

fprintf('Try with a new starting point\n')
x_0=[2;2];
[x,n]=newton(x_0,1e-6);
fprintf('The minimum of the function (x1-2)^4+(x1-2*x2)^2 is (%.2f;%.2f) \n',x(1),x(2))
fprintf('%d iterations were needed for Newtons method to converge with starting point (%.2f;%.2f)\n',n,x_0(1),x_0(1))
%'We see that in the second scenario since the first guess point is closer to the actual minimum, 
% we need one less iteration and the result is more accurate
fprintf('---------------\n')

fprintf('Exercise 2 part b) result:\n')
x_0=[2;2];
stepsize=10;
[x,n]=gradientDescent(x_0,stepsize,1e-6);
fprintf('The minimum of the function (x1-2)^4+(x1-2*x2)^2 is (%.2f;%.2f) \n',x(1),x(2))
fprintf('%d iterations were needed for Newtons method to converge with starting point (%.2f;%.2f) and %.2f stepsize \n',n,x_0(1),x_0(1),stepsize)
stepsize=0.5;
[x,n]=gradientDescent(x_0,stepsize,1e-6);
fprintf('The minimum of the function (x1-2)^4+(x1-2*x2)^2 is (%.2f;%.2f) \n',x(1),x(2))
fprintf('%d iterations were needed for gradient descent to converge with starting point (%.2f;%.2f) and %.2f stepsize\n',n,x_0(1),x_0(1),stepsize)
stepsize=0.01;
[x,n]=gradientDescent(x_0,stepsize,1e-6);
fprintf('The minimum of the function (x1-2)^4+(x1-2*x2)^2 is (%.2f;%.2f) \n',x(1),x(2))
fprintf('%d iterations were needed for gradient descent to converge with starting point (%.2f;%.2f) and %.2 fstepsize\n',n,x_0(1),x_0(1),stepsize)

%Option price function
function res=f(S,K,T,sig,r)
res=S*normcdf(d1(S,K,T,sig,r))-K.*exp(-r*T).*normcdf(d2(S,K,T,sig,r));
end

%Derivative of the option price w.r.t volatility
function res=vega(S,K,T,sig,r)
res=S*exp(-d1(S,K,T,sig,r).^2/2)./(sqrt(2*pi)).*sqrt(T);
end

function res=d1(S,K,T,sig,r)
res=(log(S/K)+(r+0.5.*sig.^2)*T)./(sig.*sqrt(T));
end

function res=d2(S,K,T,sig,r)
res=d1(S,K,T,sig,r)-sig.*sqrt(T);
end

%Computes one iteration of the Newton method
function res=iterateStep(S,K,T,sig,r)
res=sig-(f(S,K,T,sig,r)-10.78)./vega(S,K,T,sig,r);
end

function res=hessian(x)
x1=x(1);
res=[12*(x1-2)^2+2, -4; -4 ,8];
end

function res=gradient(x)
res=[ 4*(x(1)-2)^3+2*(x(1)-2*x(2)); -4*(x(1)-2*x(2)) ];
end

function [res,n_iter]=newton(x_0,eps)
x=x_0;
%Keeps track of the number of iterations needed
n=0;
%While the largest coordinate of the gradient is higher than the threshold, we keep
%iterating
while max(gradient(x))>eps
    n=n+1;
    x=x-pinv(hessian(x))*(gradient(x));
end
res=x;
n_iter=n;
end

function [res,n_iter]=gradientDescent(x_0,stepsize,eps)
x=x_0;
%Keeps track of the number of iterations needed
n=0;
%While the largest coordinate of the gradient is higher than the threshold, we keep
%iterating
while max(gradient(x))>eps
    n=n+1;
    x=x-stepsize*gradient(x);
end
res=x;
n_iter=n;
end