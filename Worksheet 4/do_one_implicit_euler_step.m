function [ T_new ] = do_one_implicit_euler_step( N_x,N_y,tau,T )
%DO_ONE_IMPLICIT_EULER_STEP computes the updated solution for the heat
%transfer equation with given mesh, timestep and initial setup using
%explicit euler method

h_x=1/(N_x+1);
h_y=1/(N_y+1);

%%T_n+1=T_n+tau*RHS(T_n+1);
warning('implementation not correct!');


RHS=zeros(N_x+2,N_y+2);

for i = 2:N_x+1
    for j = 2:N_y+1
        RHS(i,j)=(T(i-1,j)-2*T(i,j)+T(i+1,j))/h_x^2+...
            (T(i,j-1)-2*T(i,j)+T(i,j+1))/h_y^2;
    end
end

%do explicit euler step
T_new=T+tau*RHS;

end

