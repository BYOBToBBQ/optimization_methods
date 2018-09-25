S=100;
K=100;
T=1.5;
r=0.04;
eps=0.01

sigma=1;
while abs(f(S,K,T,sigma,r)-10.78)>eps
    sigma=x1(S,K,T,sigma,r);
end
fprintf('The implied volatility is %.4f \n',sigma)

function res=f(S,K,T,sig,r)
res=S*normcdf(d1(S,K,T,sig,r))-K.*exp(-r*T).*normcdf(d2(S,K,T,sig,r));
end


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
function res=x1(S,K,T,sig,r)
res=sig-(f(S,K,T,sig,r)-10.78)./vega(S,K,T,sig,r);
end





