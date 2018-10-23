
clear all;
close all;
%Define parameters
beta=0.9;
A=2;
B=1.5;
NumPoints =10;
%Discretize the state space around the optimal solution
x_bar = 0;
x_lo = x_bar-1;
x_hi = x_bar+1;
step = (x_hi-x_lo)/NumPoints;
X = x_lo:step:x_hi;
n=length(X);
%We build an n×n matrix whose columns are output at each value of X
XI = ones(n,1)*X;
%Then another n×n matrix whose columns are output at each value of X
XJ = ones(n,1)*X;
% We can construct the distance between all the Xi-Xj
X_diff=XI-XJ';
% Calculate the utility arising from each level of distance
U=-A*X_diff.^2-B*XI^2;
%Take an initial guess at the value function
V = ones(n,1);
VV=V*ones(1,n);
%Apply the operator:;
W=U+beta*VV;
V=max(W)';
%Main iteration loop for the value function.
flag=1;
while (flag > 10^(-5))
    VV=V*ones(1,n);
    W=U+beta*VV;
    V1=max(W)';
    flag = max(abs(V1-V));
    V=V1;
end
%When the value function has converged, find the policy function i.e. the Xj that gives the maximum value of the operator for each Xi. In order to accomplish this, we first find the vector of indices where W takes its maximum value:
[val,ind]=max(W);
%Then we use the indices to pick out the corresponding values of Xi.
XI_star = X(ind);

figure
plot(XI_star)
xlabel('grid')
ylabel('x')
title('Best Path to minimize the cost of travelling from -1 to 1')



