function [ T ] = initialize_Temperature_Grid(N_x, N_y)
%INITIALIZE_TEMPERATURE_GRID Summary of this function goes here
%   Detailed explanation goes here
T = ones(N_x+2,N_y+2);
T(1,:) = 0;
T(N_y+2,:) = 0;
T(:,1) = 0;
T(:,N_x+2)=0;

end

