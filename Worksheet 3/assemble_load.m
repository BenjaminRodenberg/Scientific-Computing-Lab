function [ F ] = assemble_load( N_x,N_y,f,BC)
%ASSEMBLE_LOAD Summary of this function goes here
%   Detailed explanation goes here

N=(N_x+2*BC)*(N_y+2*BC);
length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

F=zeros(N,1);

for i = 1-BC:N_x+BC
    for j = 1-BC:N_y+BC
        c_node=(i+BC)+(N_x+2*BC)*(j-1+BC);
        
        x_node=i*h_x;
        y_node=j*h_y;
        
        F(c_node)=F(c_node)+f(x_node,y_node);
    end
end

end

