function [ errors ] = calculate_analitical_errors( sols )

%get data from struct
tau_range = sols.tau;
T_end = sols.T_end;

%initialize variables
errors = zeros(1,numel(tau_range));
errors(:) = inf;

%iterate over used timestep sizes
for i = 1:numel(tau_range)
    
    %get data for current timestep size
    tau = tau_range(i);
    field_name = sprintf('tau%i',i);
    y = sols.(field_name).y;
    y_ref = sols.(field_name).y_ref;
    
    %calculate L2-Error
    curr_error = error_norm( y , y_ref , tau , T_end);
    
    %if error is equal to zero the function values are usually not defined
    %and therefore the error is inf
    if(or(curr_error == Inf,curr_error == 0))
        errors(i) = Inf;
    else
        errors(i) = curr_error;
    end
end