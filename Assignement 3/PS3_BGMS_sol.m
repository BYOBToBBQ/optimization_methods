x_0 =[2;1];
eps= 0.000001;
mu=10;
increment=1.01;

%the number of iterations is counted as the number of time we run newton's
%algorithm on a previously found minimum until convergence
[result,n_iter] =newton_penalty(mu,x_0,eps,increment);
fprintf('Increment size: %.2f\n',increment)
fprintf('Result is [%.4f,%.4f] for %d iterations\n',result(1),result(2),n_iter);
fprintf('Min point is [%.4f,%.4f]\n',result(1),result(2));
fprintf('Function evaluation at min is [%.4f]\n',fct(result));
fprintf('Equality constraint [%.4f]\n',(abs(result(2)-result(1)^2)));

%Try a faster increment
increment=10;
fprintf('We try a faster increment: %d\n',increment);
[result,n_iter] =newton_penalty(mu,x_0,eps,increment);
fprintf('Result is [%.4f,%.4f] for %d iterations\n',result(1),result(2),n_iter);
fprintf('Min point is [%.4f,%.4f]\n',result(1),result(2));
fprintf('Function evaluation at min is [%.4f]\n',fct(result));
fprintf('Equality constraint [%.4f]\n',(abs(result(2)-result(1)^2)));
%We see that a faster increment implies a faster convergence for the
%same starting mu, which is intuitive

increment=10000000000;
fprintf('We try a very large increment: %d\n',increment);
fprintf('Disclaimer, the algorithm will not converge, so feel free to stop');
[result,n_iter] =newton_penalty(mu,x_0,eps,increment);
fprintf('Result is [%.4f,%.4f] for %d iterations\n',result(1),result(2),n_iter);
fprintf('Min point is [%.4f,%.4f]\n',result(1),result(2));
fprintf('Function evaluation at min is [%.4f]\n',fct(result));
fprintf('Equality constraint [%.4f]\n',(abs(result(2)-result(1)^2)));

function [res,n_iter]= newton_penalty(mu,x_0,eps,increment)
n=0;
x_min = x_0;
%We accept a solution only if it is close to the feasible region
while (abs(x_min(2)-x_min(1)^2)>eps)
    n=n+1;
    %increase mu and run once again newton on the new penalty function
    %with as a starting guess the previously computed minimum
    mu=mu*increment;
    x_min = mu_step(mu,x_min,eps);
end
    res=x_min;
    n_iter=n;
end

%Newton step for a given mu
function res=mu_step(mu,x_0,eps)
x=x_0;
    while max(abs(gradient(mu,x)))>eps
        x=x-hessian(mu,x)\gradient(mu,x);
    end
res=x;
end

function res=fct(x)
res=(x(1)-2)^4 + (x(1)-2*x(2))^2;
end

function res=gradient(mu,x)
res=[4*(x(1)-2)^3+2*(x(1)-2*x(2))-4*mu*x(1)*(x(2)-(x(1)^2));-4*(x(1)-2*x(2))+2*mu*(x(2)-(x(1)^2))];
end

function res= hessian(mu,x)
res=[12*(x(1)-2)^2+2+8*mu*(x(1)^2)-4*mu*(x(2)-(x(1)^2)) -4-4*mu*x(1);-4-4*mu*x(1) 8+2*mu];
end