close all;
clear all;

% PDE:
% T_t = T_xx + T_yy

% Boundary Conditions
% T=0 on del[0,1]^2

% Initial Conditions
% T=1 on [0,1]^2

%--------------------------------------------------------------------------

vector_tau=[1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096];
vector_N_x=[3 7 15 31];

interesting_time=[1/8,2/8,3/8,4/8];

subplot_id=1;
FONTSIZE=12;
for j = 1:numel(interesting_time)
    figure(j)
    set(gcf,'Position',[100 100 1100 500]);
    subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
    hold on
    set(gcf,'Color','white');
    text(0.1,0,['$t = ',mat2str(interesting_time(j)),'$'],'interpreter','latex','Fontsize',FONTSIZE);
    set(gca,'Color','white');
    set(gca,'XColor','white');
    set(gca,'YColor','white');
    hold off
end

subplot_id = 1;
for i = 1:numel(vector_tau)
    subplot_id = subplot_id +1;
    for j = 1:numel(interesting_time)
        figure(j)
        subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
        hold on
        set(gcf,'Color','white');
        text(0.2,0,['$\tau = \frac{1}{',mat2str(1./vector_tau(i)),'}$'],'interpreter','latex','Fontsize',FONTSIZE);
        set(gca,'Color','white');
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        hold off
    end
end

subplot_id=1;
for i = 1:numel(vector_N_x)
    subplot_id=subplot_id+1+numel(vector_tau);
    for j = 1:numel(interesting_time)
        figure(j)
        subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
        hold on
        set(gcf,'Color','white');
        text(0.1,0.5,['$N_x = ',mat2str(vector_N_x(i)),'$'],'interpreter','latex','Fontsize',FONTSIZE);
        set(gca,'Color','white');
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        hold off
    end
end

subplot_id=length(vector_tau)+2;
%choose gridsize
for N_x=vector_N_x
    
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
    for tau=vector_tau
        
        subplot_id=subplot_id+1;
        
        % solve heat equation explicitly for given grid and timestep size
        [T,t]=solve_with_explicit_euler(N_x,N_y,tau);
        
%         % solve heat equation implicitly for given grid and timestep size
%         [T,t]=solve_with_implicit_euler(N_x,N_y,tau);
        
        
        for i=1:numel(interesting_time)
            figure(i)
            if mod(subplot_id-1,length(vector_tau)+1)==0
                subplot_id = subplot_id+1;
            end
            subplot(length(vector_N_x)+1,length(vector_tau)+1,subplot_id)
            hold on
            axis([0 1 0 1 0 1])
            %             warning('better contour or surf?');
            surf(X,Y,T(:,:,interesting_time(i)/tau+1));
            view(3);
            hold off
        end
    end
end