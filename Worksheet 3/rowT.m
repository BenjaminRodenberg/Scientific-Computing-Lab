function [ rowT ] = rowT( i, j, Nx, Ny )
%T Summary of this function goes here
%   Detailed explanation goes here

%T = sin(pi*x)*sin(pi*y);

%Ts surrounding node will be weighted:
%One row has Nx * Ny elements
%Use a vector of zeros and set the four relevant values
Tvector = zeros( Nx * Ny, 1 );
%Tset @(i,j) = Tvector((j-1)*Nx+i);

%boundaries not regarded at the moment, as they are 0 in this case

hx = 1 / ( Nx + 1 );
hy = 1 / ( Ny + 1 );

%set the 5 non-zero vector elements
%vector is adressed as: ((j-1)*Nx+i)

% Tset(i-1,j) = 1/hx^2;
% Tset(i,j) = -2*(hx^2+hy^2)/(hx^2*hy^2);
% Tset(i+1,j) = 1/hx^2;
% Tset(i,j-1) = 1/hy^2;
% Tset(i,j+1) = 1/hy^2;

i_temp_positions = [ i-1, i, i+1, i,   i ];
j_temp_positions = [   j, j,   j, j-1, j+1 ];
outcomes = [ 1/hx^2, -2*(hx^2+hy^2)/(hx^2*hy^2), 1/hx^2, 1/hy^2, 1/hy^2];

for k = 1: numel( i_temp_positions ) 
    if position_inside_boundaries( i_temp_positions(k), Nx, j_temp_positions(k), Ny )
        Tvector((j_temp_positions(k)-1)*Nx+i_temp_positions(k)) = outcomes(k);
    end
end

Tvector.'
rowT = Tvector;


end

