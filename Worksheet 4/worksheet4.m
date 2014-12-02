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

set_custom_plot(vector_N_x, vector_tau, interesting_time);

subplot_id=length(vector_tau)+2;
subplot_id_implicit=length(interesting_time)+2;

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
        
        % solve heat equation implicitly for given grid and timestep size  
        if tau == vector_tau(1);                     
            [T_implicit]=solve_with_implicit_euler(N_x,N_y,tau);
            for i=1:numel(interesting_time)
                figure(5)
                subplot_id_implicit = subplot_id_implicit+1;              
                subplot(length(vector_N_x)+1,length(interesting_time)+1,subplot_id_implicit);
                hold on
                axis([0 1 0 1 0 1])            
                surf(X,Y,T_implicit(:,:,interesting_time(i)/vector_tau(1)+1));
                view(3);
                hold off
            end   
            subplot_id_implicit = subplot_id_implicit + 1;
        end
        
        for i=1:numel(interesting_time)
            figure(i)
            if mod(subplot_id-1,length(vector_tau)+1)==0
                subplot_id = subplot_id+1;
            end
            subplot(length(vector_N_x)+1,length(vector_tau)+1,subplot_id)
            hold on
            axis([0 1 0 1 0 1])            
            surf(X,Y,T(:,:,interesting_time(i)/tau+1));
            view(3);
            hold off
        end
    end
end

for i=1:numel(interesting_time) %1. 1/8, 2. 2/8, 3. 3/8, 4. 4/8 
    figure_handle = figure(i);
    hgexport(gcf, ['ExplicitEulerForTime' num2str(i) '.jpg'],...
             hgexport('factorystyle'), 'Format', 'jpeg');    
end

