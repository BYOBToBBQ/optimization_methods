%% PROBLEM SET 01 - OPTIMIZATION METHODS
%% BRODARD LIONEL, MARCHAL ANTOINE, SCHONENBERGER ANTOINE
%% Exercice 01
% Declaration of variables
% Stock price
S=100;
% Strike price
K=100;
% Interest rate
r=0.04;
% Dividend yield
q=0.00;
% Time now
t=0;
% Maturity time
T=1.5;
% Various volatilities
sigma=(0.05:0.001:1)';
% Call prices
C=10.78;
% Implied vol
implied_sigma=impvol(C,S,K,r,t,T,q);
% Check results
error_max=max(abs(C-call(S,K,r,implied_sigma,t,T,q)));
% Declaration of functions
% Function d1
function d=d1(S,K,r,sigma,t,T)
   d=(log(S./K)+(r+sigma.^2*0.5).*(T-t))./(sigma.*sqrt(T-t));
end
% Function d2
function d=d2(S,K,r,sigma,t,T)
   d=(log(S./K)+(r-sigma.^2*0.5).*(T-t))./(sigma.*sqrt(T-t));
end
% Function normal cdf, 
function p=Phi(x)
   p=normcdf(x);
end
% Function normal pdf,
function p=PhiPrime(x)
   p=normpdf(x);
end
% Function Call Price,
function c=call(S,K,r,sigma,t,T,q)
    if nargin>6
        c=exp(-q.*(T-t)).*call(S,K,r-q,sigma,t,T);
    else
        c=S.*Phi(d1(S,K,r,sigma,t,T))-K.*exp(-r.*(T-t)).*Phi(d2(S,K,r,sigma,t,T));
    end
end
% Function Call Vega,
function v=call_vega(S,K,r,sigma,t,T,q)
    if nargin>6
        v=exp(-q.*(T-t)).*call_vega(S,K,r-q,sigma,t,T);
    else
        v=S.*PhiPrime(d1(S,K,r,sigma,t,T)).*sqrt(T-t);
    end
end

% Black-Scholes implied vol

function sigma=impvol(C,S,K,r,t,T,q,tol)

T=T-t;

if nargin<8
    tol=1e-6;
end

if nargin<7 || isempty(q)
    q=0;
end

F=S*exp((r-q).*T);
G=C.*exp(r.*T);

alpha=log(F./K)./sqrt(T);
beta=0.5*sqrt(T);

% Now we need to solve G=F Phi(d1)- K Phi(d2) where
% d1=alpha/sigma+beta and d2=alpha/sigma-beta

a=beta.*(F+K);
b=sqrt(2*pi)*(0.5*(F-K)-G);
c=alpha.*(F-K);

disc=max(0,b.^2-4*a.*c);

sigma0=(-b+sqrt(disc))./(2*a);

sigma=NewtonMethod(sigma0);

    function s1=NewtonMethod(s0)
        
        s1=s0;
        count=0;
        f=@(x) call(S,K,r,x,0,T,q)-C;
        fprime=@(x) call_vega(S,K,r,x,0,T,q);
        
        max_count=1e3;
        
        while max(abs(f(s1)))>tol && count<max_count
            count=count+1;
            
            s0=s1;
            s1=s0-f(s0)./fprime(s0);
        end
        
        if max(abs(f(s1)))>tol
            disp('Newton method did not converge')
        end
    end

end

%% Exercice 02
function f_x = fx(x_1,x_2)
f_x=(x_1-2)^4 + (x_1-2*x_2)^2;
end

