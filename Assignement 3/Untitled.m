x_0 =[3;3];
eps= 0.001;
mu=10;
increment=1;

result =newton(mu,x_0,eps,increment);
fprintf('Result is [%.4f,%.4f]',result(1),result(2));
function res= newton(mu,x_0,eps,increment)
x_last=x_0;
x_min = mu_step(mu,x_0);
if(norm(x_min-x_last)<eps)
    res=x_min;
else
    newton(mu+increment,x_0,eps,increment);
end

end
function res=mu_step(mu,x_0)
x=x_0;
    while max(abs(gradient(mu,x)))>0.00001
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