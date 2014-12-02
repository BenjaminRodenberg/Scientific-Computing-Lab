% PDE:
% T_t = T_xx + T_yy

% Boundary Conditions
% T=0 on del[0,1]^2

% Initial Conditions
% T=1 on [0,1]^2

%--------------------------------------------------------------------------
%geometry infomation

%length of the domain
length_x=1;
length_y=1;

%--------------------------------------------------------------------------
%discretization settings

%number of inner nodes
N_x=31;
N_y=31;
%resulting meshwidth
h_x=length_x/(N_x+1);
h_y=length_y/(N_y+1);
%timestep size
tau=1/64*1/2^6;
%final Time
TIME_final=4/8;

%--------------------------------------------------------------------------
%Defining relevant functions
%initialization of T (solution of PDE) and RHS (for ODE)
T=zeros(N_x+2,N_y+2);
T_store=zeros(N_x+2,N_y+2,ceil(TIME_final/tau)+1);
RHS=zeros(N_x+2,N_y+2);

%helper function returning x,y for given i,j
x_ij=@(i,j) (i-1)*h_x;
y_ij=@(i,j) (j-1)*h_x;

%Laplacian (FD)
L_T_ij=@(T,i,j) (T(i-1,j)-2*T(i,j)+T(i+1,j))/h_x^2+...
    (T(i,j-1)-2*T(i,j)+T(i,j+1))/h_y^2;

%--------------------------------------------------------------------------
%set initial conditions
for i = 2:N_x+1
    for j = 2:N_y+1
        T(i,j)=1;
    end
end

%--------------------------------------------------------------------------
%time loop
for t_current=0:tau:TIME_final
    %save data
    T_store(:,:,t_current/tau+1)=T;
    %update rhs of ODE
    for i = 2:N_x+1
        for j = 2:N_x+1
            RHS(i,j)=L_T_ij(T,i,j);
        end
    end
    %do one timestep (explicit euler)
    T=T+tau*RHS;    
end
%--------------------------------------------------------------------------
%postprocessing
X=zeros(size(T));
Y=zeros(size(T));
for i = 1:N_x+2
    for j = 1:N_y+2
        X(i,j)=x_ij(i,j);
        Y(i,j)=y_ij(i,j);
    end
end

figure(1)
hold on
view(3)
h=[];
axis([-.5 1.5 -.5 1.5 -.5 1.5])
for i=1:size(T_store,3)
    delete(h);
    h=surf(X,Y,T_store(:,:,i));  
    pause(10*tau);    
end
hold off