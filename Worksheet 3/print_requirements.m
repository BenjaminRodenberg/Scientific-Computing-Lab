function [ ] = print_requirements( method_name, N_x, solutions )

requirements = zeros( 2, numel(N_x) );

index_full_matrix   = 1;
index_sparse_matrix	= 2;
index_gauss_seidl	= 3;

for method_id = [index_full_matrix, index_sparse_matrix, index_gauss_seidl]
    fprintf('------------------------------------------------------------------------------------------------------------\n')
    fprintf('%s',(method_name{method_id}))
    
    for i=1:numel(N_x)
        current_N_x = N_x(i);
        Properties = {'Runtime';'Storage'};    
        requirements(:, i) = [solutions.(method_name{method_id}).(['N_x' num2str(current_N_x)]).runtime...
               solutions.(method_name{method_id}).(['N_x' num2str(current_N_x)]).storage]; 
    end           
    
    Nx_7 = requirements(:,1);
    Nx_15 = requirements(:,2);
    Nx_31 = requirements(:,3);
    Nx_63 = requirements(:,4);
    Requirements = table(Nx_7, Nx_15,Nx_31, Nx_63,'RowNames',Properties)
    fprintf('------------------------------------------------------------------------------------------------------------\n')
end

end
