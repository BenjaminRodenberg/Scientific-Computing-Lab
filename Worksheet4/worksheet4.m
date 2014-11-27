% clear all;
% close all;

length_x = 1;
length_y = 1;

%continuous load
f=@(x,y)-2*pi^2*sin(pi*x).*sin(pi*y);
%discrete load
with_boundaries=0;

T_analytic=@(x,y)sin(pi*x).*sin(pi*y);

solutions = struct();
%the indeces are defining the solution methods:
index_explicit_euler   = 1;
index_implicit_euler	= 2;
index_gauss_seidl	= 3;
 
method_name={'ExplicitEuler','ImplicitEuler'};

%Removed last grid size for testing purposes, as Gauss-Seidel takes forever
N_x = [3 7 15 31];
N_y = [3 7 15 31];
number_of_grid_sizes = numel(N_x);
fig_id =1;
storage = 0;
tau_range=(1/64).^(1:1:7);
sample_times = (1:4)/8;

%My interpretation of the problem: At time 0, the whole grid (boundary
%excluded) has temperature 1. Afterwards, due to the boundaries with
%temperature 0, the whole grid slowly cools down.
%Therefore, the temperature of the system converges to zero over time
%a) T(x,y,t) --t=inf--> 0
%Declining function for T

for method_id = [index_explicit_euler, index_implicit_euler, index_gauss_seidl]
%     hold on
    previous_N_x = 0;        
 
    for current_grid_index = 1:number_of_grid_sizes
        current_N_x = N_x(current_grid_index);
        current_N_y = N_y(current_grid_index);     

        T = initialize_Temperature_Grid(current_N_x, current_N_y);
        h_x = length_x/(current_N_x+1);
        h_y = length_y/(current_N_y+1);

        switch method_id
%             case index_full_matrix
%                 [T,runtime,storage] = solve_with_explicit_euler(current_N_x,...
%                                     current_N_y, b, method_name{method_id} );                
%             case index_sparse_matrix
%                 [T,runtime,storage] = solve_with_implicit_euler(current_N_x,...
%                                     current_N_y, b, method_name{method_id} );                
%             case index_gauss_seidl
%                 [T,runtime,storage] = solve_with_gauss_seidl(current_N_x,...
%                                     current_N_y, b, method_name{method_id} );
%             otherwise 
%                 disp('No solution method specified for this method_id');
%                 break
        end                  
        
    end
    
end