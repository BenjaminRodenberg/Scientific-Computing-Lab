N_x=7;
N_y=7;
N=N_x*N_y;
length_x = 1;
length_y = 1;

h_x = length_x/(N_x+1);
h_y = length_y/(N_y+1);

%continuous load
f=@(x,y)-2*pi^2*sin(pi*x).*sin(pi*y);
%discrete load
with_boundaries=0;
b=assemble_load(N_x,N_y,f,with_boundaries);

L=assemble_laplacian(N_x,N_y);


%solve System with MATLAB direct solver
T=L\b;



%corresponding grid
[x,y]=meshgrid([h_x:h_x:length_x-h_x],[h_y:h_y:length_y-h_y]);

%convert solution from vector to matrix
Z=zeros(N_y,N_x);
for i = 1:numel(x)
    T_index=get_discrete_index(x(i),y(i),h_x,h_y,N_x,N_y,with_boundaries);
    Z(i)=T(T_index);
end

T_ana=@(x,y)sin(pi*x).*sin(pi*y);
%calculate Error
E=error_norm(Z,T_ana(x,y))

%Add homogeneous boundary values
Z=[zeros(1,N_x+2);zeros(N_y,1),Z,zeros(N_y,1);zeros(1,N_x+2)];
%Create grid for plotting
[X,Y]=meshgrid([0:h_x:length_x],[0:h_y:length_y]);
[XX,YY]=meshgrid([0:.1:length_x],[0:.1:length_y]);

figure(1)
hold on
title('solution of PDE')
surf(X,Y,Z,'FaceColor','interp')
mesh(XX,YY,T_ana(XX,YY),'FaceColor','none')
xlabel('x')
ylabel('y')
hold off

figure(2)
hold on
title('solution of PDE')
contour(X,Y,Z)
xlabel('x')
ylabel('y')
hold off