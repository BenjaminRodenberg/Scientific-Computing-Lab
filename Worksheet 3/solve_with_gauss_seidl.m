function [ T, runtime, storage ] = solve_with_gauss_seidl( current_N_x, current_N_y, b, method_name )
%c) Implement a Gauss-Seidel solver for the system in a) as a function 
    %of the right hand side b, Nx, and Ny

tic;
disp(['calculating T for N= ' mat2str(current_N_x) ', method: ' method_name '...'] );
T = perform_gauss_seidl(current_N_x,current_N_y, b);
runtime = toc;
disp(['done!'] );
storage = numel(b) + numel(T); %+5 for crosshairs_stencil

end