clear all;
close all;
%Define parameters
beta=0.99;
A=8;
B=4;
%+1 to get actual number of points on grid
NumPoints =999;
%Discretize the state space around the optimal solution
x_bar = 0;
x_lo = x_bar-1;
x_hi = x_bar+1;
step = (x_hi-x_lo)/NumPoints;
X = x_lo:step:x_hi;
 
n=length(X);
%NumPoints+1, it's the actual number of points on the grid

%arrange our function as a vector since number of grid points is countable
V_0 = -(A+B)*X.^2;
M = zeros(n, n);
max_stop =1;
while max_stop > 0.00001
    V_next = zeros(1,n);
    %brute force test all possible values of the function
    for i=1:n
        for j=1:n
            M(i,j)=-A*(X(i)-X(j))^2-B*X(i)^2+beta*V_0(j);
        end
    end
    
    %select only the maximum along x(j). i.e select the max of each row to
    %get the function
    for i=1:n
        V_next(i)= max(M(i,:));
    end
    
    %check convergence condition
    max_stop = max(abs(V_next-V_0));
    V_0=V_next;
end

plot(X,V_0);
hold on
%theoretical analytically derived
a1= ((1/beta + 1 + B/A)-sqrt((1/beta + 1 + B/A)^2-4/beta))/2;

V_th = -(A*(1-a1))/(beta*a1)*X.^2;

plot(X,V_th)
hold off
%We see that both functions coincide perfectly, so much so that they are
%undistinguishable.

