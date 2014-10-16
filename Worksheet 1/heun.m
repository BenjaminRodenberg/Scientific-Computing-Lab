function [ t,f ] = heun(df,f0,tau,T_end)
% numerical solution of ODE using heun scheme

t=0:tau:T_end;

f=zeros(size(t));

f(1)=f0;

for n = 1:(numel(f)-1)
    
    f_mid=f(n)+tau*df(t(n),f(n));    
    f(n+1)=.5*(f(n)+f_mid+tau*df(t(n+1),f_mid));
    
end

end

