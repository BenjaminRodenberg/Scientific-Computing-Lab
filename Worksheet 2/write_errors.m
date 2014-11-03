function [  ] = write_errors( sols,method_name )
%WRITEERRORS writes the error and the error reduction of a numerical method
%into command window.

%get data from stuct
tau_range=sols.tau;
errors=sols.errors;
error_reduction=sols.error_reduction;
computation_time=sols.computation_time;

fprintf('%s:\n',method_name)
fprintf('------------------------------------------------------------------------------------------------------------\n')

%write used timestep size
fprintf('tau              |');
for i = 1:numel(tau_range)
    fprintf('%d\t',tau_range(i));
    if(tau_range(i)==1)
        fprintf('\t');
    end
end
fprintf('\n');

%write absolute error
fprintf('error            |');
for i = 1:numel(tau_range)
    fprintf('%d\t',errors(i));
    if(errors(i)==inf)
        fprintf('\t');
    end
end
fprintf('\n');

%write error reduction
fprintf('error reduction  |');
for i = 1:(numel(tau_range))
    fprintf('%d\t',error_reduction(i));
    if(isnan(error_reduction(i)))
        fprintf('\t');
    end
end
fprintf('\n');

%write computation time = cost
fprintf('computation_time |');
for i = 1:(numel(tau_range))
    fprintf('%d\t',computation_time(i));
    if(isnan(computation_time(i)))
        fprintf('\t');
    end
end
fprintf('\n');

fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('\n');

end

