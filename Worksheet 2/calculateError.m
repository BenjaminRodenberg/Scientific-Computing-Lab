function [ sols ] = calculateError( sols )
%CALCULATEERROR calculates the L2-Error of numerical solutions of an ODE
%using different ODE solvers. Finally also the error reduction with
%decreasing timestep size is calculated.

%get data from struct
tau_range=sols.tau;
T_end=sols.T_end;

%initialize variables
errors=zeros(1,numel(tau_range));
errors(:)=inf;

%iterate over used timestep sizes
for i = 1:numel(tau_range)
    
    %get data for current timestep size
    tau=tau_range(i);
    field_name=sprintf('tau%i',i);
    y=sols.(field_name).y;
    y_ref=sols.(field_name).y_ref;
    
    %calculate L2-Error
    curr_error = Error_norm( y , y_ref , tau , T_end);
    
    %if error is equal to zero the function values are usually not defined
    %and therefore the error is inf
    if(or(curr_error == Inf,curr_error == 0))
        errors(i) = Inf;
    else
        errors(i) = curr_error;
    end
end

%save results
sols.errors=errors;

%calculate error reduction
error_reduction=zeros(1,numel(tau_range));
error_reduction(1)=nan;

for i = 2:(numel(tau_range))
    if(errors(i)==inf||errors(i-1)==inf)
        %for inf error no reduction can be calculated
        error_reduction(i)=nan;
    else
        %reduction factor is the ration of old error and new error. If the
        %error is reduced with decreasing timestep size the reduction
        %factor is greater 1.
        error_reduction(i)=errors(i-1)/errors(i);
    end
end

%save result
sols.error_reduction=error_reduction;

end

