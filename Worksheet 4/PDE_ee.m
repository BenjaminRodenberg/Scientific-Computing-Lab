% PDE:
% T_t = T_xx + T_yy

% Boundary Conditions
% T=0 on del[0,1]^2

% Initial Conditions
% T=1 on [0,1]^2

%--------------------------------------------------------------------------

subplot_id=0;
interesting_time=[1/8,2/8,3/8,4/8];

% for fig_id = 1:4
%     figure(fig_id)
%     subplot(5,8,1)
%     title(['time = ',mat2str(interesting_time(fig_id))]);
% end

% nx_fig_id=1;
% tau_fig_id=1;

%choose gridsize
for N_x=[3 7 15 31]
    
    %     %prepare subplots
    %     nx_fig_id=nx_fig_id+8;
    %     for fig_id = 1:4
    %         figure(fig_id)
    %         subplot(5,8,nx_fig_id)
    %         title(['N_x = ',mat2str(N_x)]);
    %     end
    
    %setup grid
    N_y=N_x;
    
    X=zeros(N_x+2,N_y+2);
    Y=zeros(N_x+2,N_y+2);
    
    %helper function returning x,y for given i,j
    x_ij=@(i,j) (i-1)*1/(N_x+1);
    y_ij=@(i,j) (j-1)*1/(N_y+1);
    
    for i = 1:N_x+2
        for j = 1:N_y+2
            X(i,j)=x_ij(i,j);
            Y(i,j)=y_ij(i,j);
        end
    end
    
    %choose timestep size
    for tau=[1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096]
        
        subplot_id=subplot_id+1;
        %         %prepare subplots
        %         tau_fig_id=tau_fig_id+1;        
        %         for fig_id = 1:4
        %             figure(fig_id)
        %             subplot(5,8,tau_fig_id)
        %             title(['tau = ',mat2str(tau)]);
        %         end
        
        %solve heat equation for given grid and timestep size
        [T,t]=solve_Heat(N_x,N_y,tau);
        
        %postprocessing
        for fig_id = 1:4
            
            interesting_id = interesting_time(fig_id)/tau+1;
            
            figure(fig_id)
            subplot(4,7,subplot_id)
            hold on
            axis([0 1 0 1 0 1])
            warning('better contour or surf?');
            contour(X,Y,T(:,:,interesting_id));            
            hold off
        end
    end
end