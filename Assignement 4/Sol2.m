clear all;
close all;
%Define parameters
beta=0.95;
A=0.5;
B=7.5;
%+1 to get actual number of points on grid
NumPoints =199;
%Discretize the state space around the optimal solution
x_bar = 0;
x_lo = x_bar-1;
x_hi = x_bar+1;
step = (x_hi-x_lo)/NumPoints;
X = x_lo:step:x_hi;
 
n=length(X);
%NumPoints+1, it's the actual number of points on the grid

V_0 = -(A+B)*X.^2;
M = zeros(n, n);
max_stop =1;
while max_stop > 0.00001
    V_next = zeros(1,n);
    for i=1:n
        for j=1:n
            M(i,j)=-A*(X(i)-X(j)^2-B*X(i)^2)+beta*V_0(j);
        end
    end
    
    for i=1:n
        V_next(i)= max(M(i,:));
    end
    
    max_stop = max(abs(V_next(i)-V_0(i)));
    V_0=V_next;
end

plot(X,V_0);
