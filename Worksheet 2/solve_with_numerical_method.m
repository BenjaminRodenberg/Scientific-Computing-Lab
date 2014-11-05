function [ sol ] = solve_with_numerical_method( method_id,tau_range,T_end,p0,dp,p_ana )
%SOLVEWITHNUMERICALMETHOD solves the given ODE dp with different timesteps
%given in tau_range until the time T_end is reached. The initial value p0
%is used. The ODE solver is specified by method_ID. The reference solution
%is compared using the analytical solution p_ana of the ODE.

%initialize variables
sol=struct();

sol.tau=tau_range;

computation_time=zeros(1,numel(tau_range));

for i = 1:numel(tau_range)
    tau=tau_range(i);
        
    %chooses the method corresponding to method_id
    %for implicit methods dp has to be a symbolic function!
    switch method_id        
        case 1
            [t,y,method_computation_time] = explicit_euler(dp,p0,tau,T_end);
        case 2
            [t,y,method_computation_time] = heun(dp,p0,tau,T_end);
        case 3
            [t,y,method_computation_time] = implicit_euler(dp,p0,tau,T_end);
        case 4
            [t,y,method_computation_time] = adams_moulton(dp,p0,tau,T_end);            
        case 5
            [t,y,method_computation_time] = adams_moulton_linearisation1(dp,p0,tau,T_end);            
        case 6
            [t,y,method_computation_time] = adams_moulton_linearisation2(dp,p0,tau,T_end);            
        otherwise 
            error('unknown method id');        
    end

    computation_time(i)=method_computation_time;
    
    %save results
    field_name=sprintf('tau%i',i);
    sol.(field_name).y=y;
    sol.(field_name).y_ref=p_ana(t);
    sol.(field_name).t=t; 
    sol.(field_name).tau=tau;        
end

t_ana=0:T_end/1000:T_end;
sol.T_end=T_end;
sol.y_ana=p_ana(t_ana);
sol.t_ana=t_ana;
sol.computation_time=computation_time;

end