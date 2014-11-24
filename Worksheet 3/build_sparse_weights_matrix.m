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

sparse_i_indices = [diagonal_indices_i  i_side_indices j_side_indices identity_indices_i identity_indices_j];
sparse_j_indices = [diagonal_indices_j j_side_indices i_side_indices identity_indices_j identity_indices_i];
sparse_s_values = [diagonal_values, diagonal_side_values, diagonal_side_values, diagonal_side_values, diagonal_side_values];

% example for N_x * N_y = 3*3:
%	_ 													  _
%  |  -64    16     0    16     0     0     0     0     0  |
%  |   16   -64    16     0    16     0     0     0     0  |
%  |    0    16   -64     0     0    16     0     0     0  |
%  |   16     0     0   -64    16     0    16     0     0  |
%  |    0    16     0    16   -64    16     0    16     0  |
%  |    0     0    16     0    16   -64     0     0    16  |
%  |    0     0     0    16     0     0   -64    16     0  |
%  |    0     0     0     0    16     0    16   -64    16  |
%  |_   0     0     0     0     0    16     0    16   -64 _|
%
%Matrix is built of the block matrices:
%   
%      | B I 0 |
% 16 * | I B I |
%      | 0 I B |
%
%with each block matrix' dimension 3*3 

sparse_weights_matrix = sparse( sparse_i_indices, sparse_j_indices, sparse_s_values ,N,N,max_non_zero_elements);

return





