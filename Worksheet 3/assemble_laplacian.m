function [ L ] = assemble_laplacian( N_x,N_y )
%ASSEMBLE_LAPLACIAN Assembles the discrete Laplacian Operator.
%
%                   - DIM1 and DIM2 denote the number of discrete CELLS,
%                   - n1 and n2 denote the number of NODES inside the
%                     domain .
%                     especially NO boundary nodes! 
%                   - h the meshwidth

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

STENCIL = five_pt_stencil(h_x,h_y);

N=N_x*N_y;

L=zeros(N,N);

index_vector=[-N_x,-1,0,1,N_x];

for i = 1:N_x
    for j = 1:N_y
        c_node=i+N_x*(j-1);
        
        %decide weather node is next to boundary node
        x_up=i<N_x;
        x_dwn=i>1;
        y_up=j<N_y;
        y_dwn=j>1;
        
        c_index_vector=c_node+index_vector(logical([y_dwn,x_dwn,1,x_up,y_up]));
        
        %assemble L                                 
        L(c_node,c_index_vector)=L(c_node,c_index_vector)+STENCIL(logical([y_dwn,x_dwn,1,x_up,y_up]));
    end
end

return

function [ STENCIL ] = five_pt_stencil( h_x,h_y )
%FIVE_PT_STENCIL Summary of this function goes here
%   Detailed explanation goes here

STENCIL=1/(h_x^2*h_y^2)*[h_x^2 h_y^2 -2*(h_y^2+h_y^2) h_y^2 h_x^2];

return





