function [ t,y,computation_time ] = implicit_calculation( G, dG, y0, tau, T_end , calling_function )
%Calculation for any implicit funcion whose G and dG are already kwown

t = 0:tau:T_end;
y = zeros( size(t) );
y(1) = y0;
desired_accuracy = 10^-4;

tic;
for n = 1:( numel(y) - 1 )  
    %calculate root equation
    current_G = @(x)G(x,y(n));
    current_dG = @(x)dG(x,y(n));
    
    [y(n+1),x_history] = newton( current_G, current_dG, y(n), desired_accuracy , [calling_function,' with tau = ',mat2str(tau), ' at t = ',mat2str(t(n))]);    
    
end
computation_time = toc;
end