function [ t,y,computation_time ] = heun(f,y0,tau,T_end)
% numerical solution of ODE using heun scheme

tic;

t = 0:tau:T_end;
y = zeros(size(t));
y(1) = y0;

tic;
for n = 1:(numel(y)-1)
    
    y_mid = y(n) + tau*f( y(n) );    
    y(n+1) = .5*( y(n) + y_mid + tau*f(y_mid) );
    
end
computation_time=toc;

computation_time = toc;

end

