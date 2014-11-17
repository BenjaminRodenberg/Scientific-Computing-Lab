function [ sparse_weights_matrix ] = build_sparse_weights_matrix( N_x,N_y )
% Assembles the discrete Laplacian Operator.
%
%                   - DIM1 and DIM2 denote the number of discrete CELLS,
%                   - n1 and n2 denote the number of NODES inside the
%                     domain .
%                     especially NO boundary nodes! 
%                   - h the meshwidth

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

crosshairs_stencil = create_crosshairs_stencil(h_x,h_y);

N = N_x * N_y;

max_non_zero_elements = N_x*N_y*5;
sparse_weights_matrix = sparse([],[],[],N,N,max_non_zero_elements);
%Yet to be optimized: Remove indexing of sparse matrix in the loop

index_displacement=[-N_x,-1,0,1,N_x];

for i = 1:N_x
    for j = 1:N_y
        current_node_index=i+N_x*(j-1);
        
        %decide weather node is next to boundary node
        x_up = i<N_x;
        x_dwn = i>1;
        y_up = j<N_y;
        y_dwn = j>1;
        
        neighbors_indices = current_node_index + index_displacement(logical([y_dwn,x_dwn,1,x_up,y_up]));
        
        %assemble weights_matrix         
        sparse_weights_matrix(current_node_index,neighbors_indices) = sparse_weights_matrix(current_node_index,neighbors_indices)...
                            +crosshairs_stencil(logical([y_dwn,x_dwn,1,x_up,y_up]));
    end
end

return





