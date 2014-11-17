function [ T, runtime, storage ] = solve_with_sparse_matrix( current_N_x, current_N_y, b, method_name )

tic;
sparse_weights_matrix = build_sparse_weights_matrix(current_N_x,current_N_y);
disp(['calculating T for N= ' mat2str(current_N_x) ', method: ' method_name]);
T = sparse_weights_matrix \ b;
runtime = toc;
storage = 0;%nzn(sparse_weights_matrix) + numel(T) + numel(b);
disp('done!');

end