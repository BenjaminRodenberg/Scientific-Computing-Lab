function [ sparse_weights_matrix ] = build_sparse_weights_matrix( N_x,N_y )
%Build sparse weights matrix

N = N_x * N_y;

max_non_zero_elements = N_x*N_y*5;

diagonal_indices_i = (1:N);
diagonal_indices_j = (1:N);
diagonal_values = ones(1, N) * ((N_x+1)*(N_y+1))*(-4);

[i_side_indices, j_side_indices] = get_side_values_indices( N_x, N_y );
diagonal_side_values = ones(1, N_x*(N_y-1))*((N_x+1)*(N_y+1));

identity_indices_i = (1: N_x*(N_y-1));
identity_indices_j = (N_x+1: N_x*N_y);

sparse_weights_matrix = sparse([diagonal_indices_i  i_side_indices j_side_indices identity_indices_i identity_indices_j],...    
    [diagonal_indices_j j_side_indices i_side_indices identity_indices_j identity_indices_i], ...
    [diagonal_values diagonal_side_values diagonal_side_values diagonal_side_values diagonal_side_values] ...
    ,N,N,max_non_zero_elements);

return





