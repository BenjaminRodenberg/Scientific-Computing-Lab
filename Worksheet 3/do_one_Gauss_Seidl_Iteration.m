function [ x ] = do_one_Gauss_Seidl_Iteration( N_x,N_y,b,x,BC )
%DO_GAUSS_SEIDL_ITERATION Summary of this function goes here
%   Detailed explanation goes here

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

STENCIL = five_pt_stencil(h_x,h_y);

N=(N_x+2*BC)*(N_y+2*BC);

L=zeros(N,N);

%saves displacement to neighbouring nodes
displacement_vector=[-N_x,-1,0,1,N_x];

for i = 1:N_x+2*BC
    for j = 1:N_y+2*BC
        current_node_index=i+(N_x+2*BC)*(j-1);
        
        %decide weather node is next to boundary node
        x_up=i<N_x+2*BC;
        x_dwn=i>1;
        y_up=j<N_y+2*BC;
        y_dwn=j>1;
        
        decider=logical([y_dwn,x_dwn,0,x_up,y_up]);
        current_index_vector=current_node_index+displacement_vector(decider);
        
        %local residual set to zero: 0 = b-weights*neighbouring Temperature
        %This can be solved for the temperature on the current node.
        x(current_node_index)=(b(current_node_index)-STENCIL(decider)*x(current_index_vector))/STENCIL((1+end)/2);       
    end
end

return

function [ STENCIL ] = five_pt_stencil( h_x,h_y )
%FIVE_PT_STENCIL Summary of this function goes here
%   Detailed explanation goes here

STENCIL=1/(h_x^2*h_y^2)*[h_x^2 h_y^2 -2*(h_y^2+h_y^2) h_y^2 h_x^2];

return

