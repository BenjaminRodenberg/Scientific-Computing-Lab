function [ index ] = get_discrete_index( x,y,h_x,h_y,N_x,N_y,BC )
%GET_T_VALUE Summary of this function goes here
%   Detailed explanation goes here

i_x=round(x/h_x)+BC;
i_y=round(y/h_y)+BC;
index = i_x+(i_y-1)*(N_x+2*BC);


end

