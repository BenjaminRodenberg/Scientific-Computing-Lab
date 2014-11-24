function [ T, storage ] = perform_gauss_seidl(N_x,N_y, b)
%close all;

tolerance=10^-4;
with_boundaries=0;
residual=inf;

storage = 0;

N=(N_x+2*with_boundaries)*(N_y+2*with_boundaries);

length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

T=zeros(N,1);

%continuous load
f=@(x,y)-2*pi^2*sin(pi*x).*sin(pi*y);
%discrete load
b=build_solution_vector(N_x,N_y,f,with_boundaries);

while(residual>tolerance)
    [T, residual]=do_one_Gauss_Seidl_Iteration(N_x,N_y,b,T,with_boundaries,h_x,h_y);
end