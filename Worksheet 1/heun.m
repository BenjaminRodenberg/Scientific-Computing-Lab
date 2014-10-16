function [ t,y ] = heun(f,y0,tau,T_end)
% numerical solution of ODE using heun scheme

t=0:tau:T_end;

y=zeros(size(t));

y(1)=y0;

for n = 1:(numel(y)-1)
    
    y_mid=y(n)+tau*f(t(n),y(n));    
    y(n+1)=.5*(y(n)+y_mid+tau*f(t(n+1),y_mid));
    
end

end

