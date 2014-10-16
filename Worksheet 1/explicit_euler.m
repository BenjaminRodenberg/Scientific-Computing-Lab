function [ t,f ] = explicit_euler(df,f0,tau,T_end)
% numerical solution of ODE using explicit euler scheme
t=0:tau:T_end;
f=zeros(size(t));
f(1)=f0;
for n = 1:(numel(f)-1)    
    f(n+1)=f(n)+tau*df(t(n),f(n));    
end
end