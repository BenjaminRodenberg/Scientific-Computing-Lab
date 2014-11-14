function [ answer ] = position_inside_boundaries( i_temp, Nx, j_temp, Ny )

answer = ( i_temp > 0 ) & ( i_temp < Nx+1 ) & ( j_temp>0 ) & ( j_temp < Ny+1 )

end