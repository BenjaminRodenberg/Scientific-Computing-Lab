function [ t,y ] = implicit_euler( sym_f,y0,tau,T_end )
% numerical solution of ODE using implicit euler scheme

[G,dG] = derive_G_implicit_euler(sym_f,tau);
[t, y, computation_time] = implicit_calculation( G, dG, y0, tau, T_end );
end

