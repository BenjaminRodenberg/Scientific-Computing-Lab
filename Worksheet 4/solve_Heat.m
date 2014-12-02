function [ T_store,t_store ] = solve_Heat( N_x,N_y,tau )
%SOLVE_HEAT Summary of this function goes here
%   Detailed explanation goes here

%final Time
TIME_final=4/8;

%--------------------------------------------------------------------------

%initialization of T (solution of PDE) and RHS (for ODE)
T=zeros(N_x+2,N_y+2);
T_store=zeros(N_x+2,N_y+2,ceil(TIME_final/tau)+1);

%initialization of Matrices for storing results
t_store=zeros(ceil(TIME_final/tau)+1);

%--------------------------------------------------------------------------

%set initial conditions
for i = 2:N_x+1
    for j = 2:N_y+1
        T(i,j)=1;
    end
end

%--------------------------------------------------------------------------
%time loop
for t_current=0:tau:TIME_final
    %save data
    T_store(:,:,t_current/tau+1)=T;
    t_store(t_current/tau+1)=t_current;
    %update rhs of ODE
    T = do_one_explicit_euler_step( N_x,N_y,tau,T );
end

end

