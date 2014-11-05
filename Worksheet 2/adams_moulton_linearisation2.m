function [ t,y,computation_time ] = adams_moulton_linearisation2( sym_f,y0,tau,T_end )
% numerical solution of ODE using implicit adams moulton in a linearized
% form

[G,dG] = derive_G_implicit_adams_moulton_linearisation2(sym_f,tau);

[t, y, computation_time] = implicit_calculation( G, dG, y0, tau, T_end, 'adams-moulton-lin2' );

end