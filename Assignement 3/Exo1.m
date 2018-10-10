x_0 =[20000;-20000];
eps= 0.000000001;
mu=10;
increment=1;

[result,niter] =newton(mu,x_0,eps,increment);
fprintf('Result is [%.4f,%.4f] for %d iterations',result(1),result(2),niter);
function [res,n_iter]= newton(mu,x_0,eps,increment)
n=1
x_last=x_0;
x_min = mu_step(mu,x_0);
fprintf('Min is [%.6f,%.6f]\n',x_min(1),x_min(2))
%probl�me convergence
while((max(abs(x_min-x_last)))>eps)
%while(n<100000)
    n=n+1;
    x_last=x_min;
    x_min = mu_step(mu+increment,x_min);
    fprintf('Min is [%.6f,%.6f]\n',x_min(1),x_min(2));
    fprintf('Min is %.6f\n',(max(abs(x_min-x_last))));
end
    res=x_min;
    n_iter=n;
end
function res=mu_step(mu,x_0)
x=x_0;
    while max(abs(gradient(mu,x)))>0.0000001
        x=x-pinv(hessian(mu,x))*(gradient(mu,x));
    end
res=x;
end

function res=gradient(mu,x)
res=[4*(x(1)-2)^3+2*(x(1)-2*x(2))-4*mu*x(1)*(x(2)-x(1)^2);-4*(x(1)-2*x(2))+2*mu*(x(2)-x(1)^2)];
end

function res= hessian(mu,x)
res=[12*(x(1)-2)^2+2*x(1)+8*mu*x(1)^2-4*mu*(x(2)-x(1)^2),-4-4*mu*x(1);-4-4*mu*x(1),8+2*mu];
end