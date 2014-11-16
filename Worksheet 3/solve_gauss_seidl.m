function [ T, storage ] = solve_gauss_seidl(N_x,N_y, b)


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

surface_handle=[];

iteration_counter=0;

while(residual>tolerance)
    T=do_one_Gauss_Seidl_Iteration(N_x,N_y,b,T,with_boundaries);
    residual=calculate_residual(N_x,N_y,b,T,with_boundaries);
    
%     %corresponding grid
%      [grid_x,y]=meshgrid([~BC*h_x:h_x:length_x-~BC*h_x],[~BC*h_y:h_y:length_y-~BC*h_y]);
%     
%     %convert solution from vector to matrix
%     Z=zeros(N_y+2*BC,N_x+2*BC);
% 
%     for i = 1:numel(grid_x)
%         T_index=get_discrete_index(grid_x(i),y(i),h_x,h_y,N_x,N_y,BC);
%         Z(i)=T(T_index);
%     end
    
%     T_analytic=@(x,y)sin(pi*x).*sin(pi*y); 
%     
%     if(~BC)
%         %Add homogeneous boundary values
%         Z=[zeros(1,N_x+2);zeros(N_y,1),Z,zeros(N_y,1);zeros(1,N_x+2)];
%     end
%     %Create grid for plotting
%     [X,Y]=meshgrid([0:h_x:length_x],[0:h_y:length_y]);
%     [XX,YY]=meshgrid([0:.1:length_x],[0:.1:length_y]);
%     
%     iteration_counter=iteration_counter+1;
%     
%     figure(1)
%     hold on
%     title('solution of PDE')
%     
%     delete(surface_handle);
%     surface_handle=surf(X,Y,Z,'FaceColor','interp');
%     
%     view(3);
%     mesh(XX,YY,T_analytic(XX,YY),'FaceColor','none')
%     xlabel('x')
%     ylabel('y')
%     
%     pause(.01);
%     hold off
        
    % figure(2)
    % hold on
    % title('solution of PDE')
    % contour(X,Y,Z)
    % xlabel('x')
    % ylabel('y')
    % hold off
end

% disp(['It took ',mat2str(iteration_counter),' iterations to reach an accuracy level of ',mat2str(tolerance),'.']);