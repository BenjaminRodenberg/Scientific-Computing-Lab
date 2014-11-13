%Worksheet 3

%2D stationary heat equation 

%syms sym_Txx sym_Tyy;
%syms x y;
%sym_static_equation(x,y) = -2*pi.^2*sin(pi*x)*sin(pi*y);
%static_equation = matlabFunction( simplify(  sym_stat_equation));

Ny=3;
Nx=3;

[matrix_A, matrix_x, matrix_b] = create_Matrices(Nx,Ny)

