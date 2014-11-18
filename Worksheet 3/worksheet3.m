clear all;
close all;

length_x = 1;
length_y = 1;

%continuous load
f=@(x,y)-2*pi^2*sin(pi*x).*sin(pi*y);
%discrete load
with_boundaries=0;

T_analytic=@(x,y)sin(pi*x).*sin(pi*y);

solutions = struct();
%the indeces are defining the solution methods:
index_full_matrix   = 1;
index_sparse_matrix	= 2;
index_gauss_seidl	= 3;

method_name={'FullMatrix','SparseMatrix','GaussSeidl'};

%Removed last grid size for testing purposes, as Gauss-Seidel takes forever
N_x = [7 15 31 63]; %127];
N_y = [7 15 31 63]; %127];
number_of_grid_sizes = numel(N_x);
fig_id =1;
storage = 0;

for method_id = [index_full_matrix, index_sparse_matrix, index_gauss_seidl]
    hold on
    for current_grid_index = 1:number_of_grid_sizes
        current_N_x = N_x(current_grid_index);
        current_N_y = N_y(current_grid_index);
        %b) Implement a function creating the matrix from a) as a function of Nx
        %and Ny
        b=build_solution_vector(current_N_x,current_N_y,f,with_boundaries);
        h_x = length_x/(current_N_x+1);
        h_y = length_y/(current_N_y+1);
%       d) Solve the system from a)
%       1) storing the system matrix as a normal (full) N N matrix.
%       2) storing the system matrix as a sparse matrix.
%       3) without storing the system matrix (use Gauss-Seidel).
        switch method_id
            case index_full_matrix
                [T,runtime,storage] = solve_with_full_matrix(current_N_x,...
                                    current_N_y, b, method_name{method_id} );                
            case index_sparse_matrix
                [T,runtime,storage] = solve_with_sparse_matrix(current_N_x,...
                                    current_N_y, b, method_name{method_id} );                
            case index_gauss_seidl
                [T,runtime,storage] = solve_with_gauss_seidl(current_N_x,...
                                    current_N_y, b, method_name{method_id} );
            otherwise 
                disp('No solution method specified for this method_id');
                break
        end                  
        
        %corresponding grid        
        [x,y]=meshgrid(h_x:h_x:length_x-h_x,h_y:h_y:length_y-h_y);
               
        %convert solution from vector to matrix
        Z=zeros(current_N_y,current_N_x);
        for i = 1:numel(x)
            T_index=get_discrete_index(x(i),y(i),h_x,h_y,current_N_x,current_N_y,with_boundaries);
            Z(i)=T(T_index);
        end        
        
        %calculate Error
        E=error_norm(Z,T_analytic(x,y));        
        solutions.(method_name{method_id}).(['N_x' num2str(current_N_x)]).runtime = runtime ;
        solutions.(method_name{method_id}).(['N_x' num2str(current_N_x)]).storage = storage;
        solutions.(method_name{method_id}).(['N_x' num2str(current_N_x)]).error = E;
        
        %Create grid for plotting
        [X_grid,Y_grid]=meshgrid(0:h_x:length_x,0:h_y:length_y);
        %[XX,YY]=meshgrid(0:.1:length_x,0:.1:length_y);

        if current_grid_index < 5 %no visualization for Nx = Ny = 127
            %e) Visualize the solutions as a
            %1) a coloured surface
            %2) a contour plot.
            plot_results_for_method(Z, X_grid, Y_grid, fig_id, current_grid_index,...
                                    current_N_x, current_N_y, method_name{method_id});
        end
                    
    end
    hold off 
    fig_id = fig_id + 1;
    
    %To do:
%     -Maybe find a way to allocate/build sparse matrix
%     -Structure code into subsections a)...f)

end

%f) Compare the runtimes and the storage requirements (measured by the number of
% entries of the arrays and/or vectors needed) 
print_requirements( method_name, N_x, solutions )

