function [ t,y ] = runge_kutta_4( f,y0,tau,T_end )
% numerical solution of ODE using runge-kutta scheme of fourth order

t=0:tau:T_end;

y=zeros(size(t));

y(1)=y0;

for n = 1:(numel(y)-1)
    
    dy1=f(t(n),y(n));
    dy2=f(t(n)+tau/2,y(n)+tau/2*dy1);
    dy3=f(t(n)+tau/2,y(n)+tau/2*dy2);
    dy4=f(t(n)+tau,y(n)+tau*dy3);
    
    y(n+1)=y(n)+tau*1/6*(dy1+2*dy2+2*dy3+dy4);
    
end

end

