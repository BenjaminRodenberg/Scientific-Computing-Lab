function [ ] = plot_results_for_method( Z, X_grid, Y_grid, fig_id, current_grid_index,...
                                        current_N_x, current_N_y, method_name )

%Add homogeneous boundary values
Z=[zeros(1,current_N_x+2);zeros(current_N_y,1),Z,zeros(current_N_y,1);zeros(1,current_N_x+2)];

%open figure
figure(fig_id);                
set(fig_id, 'Position', [100 100 1100 500])

subplot(2,4,(2*current_grid_index-1));
surf(X_grid,Y_grid,Z,'FaceColor','interp')
%mesh(XX,YY,T_analytic(XX,YY),'FaceColor','none')
str_title = ['                   ' method_name '- Surface and Contour: N_x, Ny = '...
            mat2str(current_N_x) ];
title(str_title);        
xlabel('x')
ylabel('y')               

subplot(2,4,(2*current_grid_index));
contour(X_grid,Y_grid,Z)        
%str_title = [method_name{method_id} '-Contour: N_x, Ny = ' mat2str(current_N_x) ];
%title(str_title);
xlabel('x')
ylabel('y')    
        
end

