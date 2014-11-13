function [matrix_A, matrix_x, matrix_b ] = create_Matrices( Nx,Ny )
%CREATE_MATRICES Summary of this function goes here
%   Detailed explanation goes here

matrix_A = zeros(Nx*Ny); 
for j = 1:1:Ny
    for i = 1:1:Nx
        vectorRow =rowT(i,j,Nx,Ny);
        matrix_A(((j-1)*Ny+i),:)= vectorRow;
    end
end



matrix_b = zeros(Nx*Ny,1);
for j = 1:1:Ny
    for i = 1:1:Nx
        x=i/(Nx+1);
        y=j/(Ny+1);
        matrix_b(((j-1)*Ny+i),1)= -2*pi*sin(pi*x)*sin(pi*y);
    end
end


matrix_x = ones(Nx*Ny,1);

end

