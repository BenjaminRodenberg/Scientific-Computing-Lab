function [ t,y ] = implicit_euler( sym_f,y0,tau,T_end )
% numerical solution of ODE using implicit euler scheme
t=0:tau:T_end;
y=zeros(size(t));
y(1)=y0;

%derive root equation for newton method symbolically
[G,dG]=derive_G_implicit_euler(sym_f,tau);
%G(y_{n+1},y_{n})=...
%dG/dy_{n+1}(y_{n+1},y_{n})=...

for n = 1:(numel(y)-1)  
    %calculate root equation
    current_G=@(x)G(x,y(n));
    current_dG=@(x)dG(x,y(n));
    %calculate root using newtons method
    y(n+1)=newton(current_G,current_dG,y(n),10^-4);    
end
end

