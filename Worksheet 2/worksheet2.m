clear all;
close all;

%symbolic formulation of the ODE dp(p)=7*(1-p/10)*p
syms sym_p;
sym_dp(sym_p)=7*(1-sym_p/10)*sym_p;


%conversion from symbolic form to matlab function form:
%i.e.: dp=@(p)7*(1-p/10).*p;
dp=matlabFunction(sym_dp);
p0=20;

T_end=5;

%analytical solution of ODE
p_ana=@(t)200./(20-10*exp(-7*t));
t_ana=0:.001:T_end;

%ODE from worksheet 1
sym_dp(sym_p)=(1-sym_p/10).*sym_p;
p0=1;
p_ana=@(t)10./(1+9*exp(-t));
dp=matlabFunction(sym_dp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) Plot the function p(t) in a graph

figure(1)
hold on
title('exact solution of ODE dp(p)=7*(1-p/10)*p')
plot(t_ana,p_ana(t_ana));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the indeces are defining the following numerical methods:
explicit_euler_index    = 1;
heun_index              = 2;
implicit_euler_index    = 3;
adams_moulton_index     = 4;
adams_moulton_l1_index  = 5;
adams_moulton_l2_index  = 6;

method_name={'explicit_euler','heun','implicit_euler',...
    'adams_moulton','adams_moulton_lin1','adams_moulton_lin2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b) Reuse the Euler method and the method of Heun ...

%different timestep size
tau_range= (.5).^(0:1:5);

% only explicit euler method and heun
for i = 1:2
    
    numerical_solutions.(method_name{i})=solveWithNumericalMethod(...
        i,tau_range,T_end,p0,dp,p_ana);         

    plotNumericalSolutions( numerical_solutions,method_name{i},i+1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d) For both methodes implemented ...

%different timestep size
clear tau_range;
tau_range=(.5).^(1:1:5);

% only implicit euler method
for i = 3
    
    numerical_solutions.(method_name{i})=solveWithNumericalMethod(...
        i,tau_range,T_end,p0,sym_dp,p_ana);         

    plotNumericalSolutions( numerical_solutions,method_name{i},i+1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g) Compare the results of the implicit methods...

%calculate error and error reduction and write results down
for i = 1:3
    numerical_solutions.(method_name{i})=...
        calculateError(numerical_solutions.(method_name{i}));
    
    writeErrors(numerical_solutions.(method_name{i}),method_name{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%