function [ T, runtime, storage ] = solve_with_full_matrix( current_N_x, current_N_y, b, method_name )

tic;
weights_matrix = build_weights_matrix(current_N_x,current_N_y);
disp(['calculating T for N= ' mat2str(current_N_x) ', method: ' method_name]);
T = weights_matrix \ b;
runtime = toc;
storage = numel(weights_matrix) + numel(T) + numel(b);
disp('done!');

end

