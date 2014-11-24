function [ crosshairs_stencil ] = create_crosshairs_stencil( h_x,h_y )
%FIVE_PT_STENCIL Summary of this function goes here
%   Detailed explanation goes here

crosshairs_stencil=(1/(h_x^2*h_y^2))*[h_x^2, h_y^2, -2*(h_y^2+h_x^2), h_y^2, h_x^2];

return
