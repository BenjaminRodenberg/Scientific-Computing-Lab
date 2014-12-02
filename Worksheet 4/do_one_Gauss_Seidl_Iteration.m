function [ A ] = do_one_Gauss_Seidl_Iteration( N_x,N_y,A,B,STENCIL )
%DO_GAUSS_SEIDL_ITERATION Summary of this function goes here
%   Detailed explanation goes here

warning('implementation not correct!');

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

%saves displacement to neighbouring nodes

for i = 2:N_x+1
    for j = 2:N_y+1
        %local residual set to zero: 0 = b-weights*neighbouring Temperature
        %This can be solved for the temperature on the current node.
        A(i,j)=(B(i,j)-(...
            A(i+1,j)*STENCIL(2,3)+...
            A(i-1,j)*STENCIL(2,1)+...
            A(i,j+1)*STENCIL(1,2)+...
            A(i,j-1)*STENCIL(3,2)))...
            /(STENCIL(2,2));        
    end
end

return

