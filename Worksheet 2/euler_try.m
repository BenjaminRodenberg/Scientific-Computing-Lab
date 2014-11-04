function [ t,y,computation_time ] = euler_try( sym_f,y0,tau,T_end )
t=0:tau:T_end;
y=zeros(size(t));
y(1)=y0;

for n = 1:(numel(y)-1)    
    y(n+1)=y(n)+tau*sym_f(t);
    y(n+1)=y(n)+tau*sym_f(y(n+1));
end


end