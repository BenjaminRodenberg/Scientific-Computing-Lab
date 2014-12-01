function [ T ] = explicit_euler( Nx, Ny, tau, T )
%explicit_euler
%Nx and Ny: Grid sizes
%Tau: delta-time
%Initial T

T_old = T;

hx = 1/(Nx + 1);
hy = 1/(Ny + 1);

weights = [1/hx^2, 1/hy^2, -2*(hx^2+hy^2)/(hx^2*hy^2)];

for i=2:Nx+1
    for j=2:Ny+1        
        T(i,j) = T_old(i,j) + tau*(...
        weights(1)*(T_old(i-1,j)+ T_old(i+1,j)) +...
        weights(2)*(T_old(i,j-1)+ T_old(i,j+1)) +...
        weights(3)*T_old(i,j));
    end
end

surf((0:hx:1), (0:hy:1), T);

end