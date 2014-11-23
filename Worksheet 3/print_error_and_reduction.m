function [ ] = print_error_and_reduction( method_name, N_x, solutions )
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('%s',(method_name))
   
Errors = zeros(2, numel(N_x));

for i=1:numel(N_x)
    current_N_x = N_x(i);
    Properties = {'Error';'Error_Reduction'};        
    Errors(:, i) = [solutions.(method_name).(['N_x' num2str(current_N_x)]).error...
           solutions.(method_name).(['N_x' num2str(current_N_x)]).error_reduction]; 
end           
    Nx_7 = Errors(:,1);    
    Nx_15 = Errors(:,2);
    Nx_31 = Errors(:,3);  
    Nx_63 = Errors(:,4);
    %Nx_127 = Errors(:,5);
    Error = table(Nx_7, Nx_15, Nx_31, Nx_63,... %Nx_127,...
    'RowNames',Properties)
    fprintf('------------------------------------------------------------------------------------------------------------\n')
end
