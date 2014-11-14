function [ T ] = do_one_Gauss_Seidl_Iteration( N_x,N_y,b,T,BC )
%DO_GAUSS_SEIDL_ITERATION Summary of this function goes here
%   Detailed explanation goes here

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

crosshairs_stencil = create_crosshairs_stencil(h_x,h_y);

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
        
        mask=logical([y_dwn,x_dwn,0,x_up,y_up]);
        masked_index_vector=current_node_index+displacement_vector(mask);
        masked_crosshairs_stencil = crosshairs_stencil(mask);
        
        %local residual set to zero: 0 = b-weights*neighbouring Temperature
        %This can be solved for the temperature on the current node.
        center_of_stencil = (numel(crosshairs_stencil)+1)/2;
        T(current_node_index)=(b(current_node_index)-masked_crosshairs_stencil*T(masked_index_vector))/crosshairs_stencil(center_of_stencil);       
    end
end

return

