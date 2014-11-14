function [ b ] = build_solution_vector( N_x,N_y,f,BC)
%ASSEMBLE_LOAD Summary of this function goes here
%   Detailed explanation goes here

N=(N_x+2*BC)*(N_y+2*BC);
length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

b=zeros(N,1);

for i = 1-BC:N_x+BC
    for j = 1-BC:N_y+BC
        current_node=(i+BC)+(N_x+2*BC)*(j-1+BC);
        
        x_node_coordinate=i*h_x;
        y_node_coordinate=j*h_y;
        
        b(current_node)=b(current_node)+f(x_node_coordinate,y_node_coordinate);
    end
end

end

