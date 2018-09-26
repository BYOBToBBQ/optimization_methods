[X,Y] = meshgrid(-2:.2:2); 
Z = (X-2)^4+(X-2*Y)^2;
surf(X,Y,Z)


x=[3;3];
n=0;
while abs(x(1))>1e-3 &&  abs(x(2))>1e-3
    n=n+1;
    x=iterate(x);
end

function res=iterate(x)
res=x-pinv(hessian(x))*(gradient(x));
end

function res=f(x1,x2)
res=(x1-2)^4+(x1-2*x2)^2;
end

function res=hessian(x)
x1=x(1);
res=[12*(x1-2)^2+2, -4; -4 ,8];
end

function res=gradient(x)
res=[ 4*(x(1)-2)^3+2*(x(1)-2*x(2)); -4*(x(1)-2*x(2)) ];
end