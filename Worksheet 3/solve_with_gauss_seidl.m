function [ T, runtime, storage ] = solve_with_gauss_seidl( current_N_x, current_N_y, b, method_name )

tic;
disp(['calculating T for N= ' mat2str(current_N_x) ', method: ' method_name '...'] );
T = perform_gauss_seidl(current_N_x,current_N_y, b);
runtime = toc;
storage = numel(b) + numel(T); %+5 for crosshairs_stencil

end