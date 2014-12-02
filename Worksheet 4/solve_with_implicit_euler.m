function [ T_history,time_history ] = solve_with_implicit_euler( N_x,N_y,tau  )
%SOLVE_WITH_IMPLICIT_EULER Summary of this function goes here
%   Detailed explanation goes here

%final Time
TIME_FINAL=4/8;

%--------------------------------------------------------------------------

%initialization of T (solution of PDE) and RHS (for ODE)
T=zeros(N_x+2,N_y+2);
T_history=zeros(N_x+2,N_y+2,ceil(TIME_FINAL/tau)+1);

%initialization of Matrices for storing results
time_history=zeros(ceil(TIME_FINAL/tau)+1);

%--------------------------------------------------------------------------

%set initial conditions, homogeneous Dirichlet boundary conditions
for i = 2:N_x+1
    for j = 2:N_y+1
        T(i,j)=1;
    end
end


tolerance = 1e-6;

%--------------------------------------------------------------------------
%time loop
for t_current=0:tau:TIME_FINAL-tau
    %save data
    T_history(:,:,t_current/tau+1)=T;
    time_history(t_current/tau+1)=t_current;
    
    %update rhs of ODE
    residual = inf;
    T_old = T;
    initial_guess = T;
    T = initial_guess;
    tic;
    while (residual > tolerance)
        [T, residual] = do_one_implicit_euler_step( N_x,N_y,tau,T, T_old );       
    end
    time = toc;
    %disp(['Time at tau=', mat2str(tau), ': ' mat2str(time)]);  

    %save final data
    T_history(:,:,TIME_FINAL/tau+1)=T;
    time_history(TIME_FINAL/tau+1)=t_current;

end


end

