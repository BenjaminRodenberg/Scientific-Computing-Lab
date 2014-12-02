function [ T_new, residual ] = do_one_implicit_euler_step( N_x,N_y,tau,T_new,T_old )
%DO_ONE_IMPLICIT_EULER_STEP computes the updated solution for the heat
%transfer equation with given mesh, timestep and initial setup using
%explicit euler method

h_x=1/(N_x+1);
h_y=1/(N_y+1);

h_x_pow2 = h_x*h_x;
h_y_pow2 = h_y*h_y;
prod_h_xy_pow2 = h_x_pow2*h_y_pow2;
global_factor = 1/(prod_h_xy_pow2+2*tau*(h_x_pow2+h_y_pow2));


for i = 2:N_x+1
    for j = 2:N_y+1
        T_new(i,j)=global_factor*(T_old(i,j)*prod_h_xy_pow2+tau*(h_y_pow2*(T_new(i-1,j)+T_new(i+1,j))+h_x_pow2*(T_new(i,j-1)+T_new(i,j+1))));
    end
end
% T_new=global_factor*T_new;

%Residual
residual = 0;
for i = 2:N_x+1
    for j = 2:N_y+1
        residual = residual + ...
            (T_new(i,j) - global_factor*...
                (T_old(i,j)*prod_h_xy_pow2+...
                    tau*(h_y_pow2*(T_new(i-1,j)+T_new(i+1,j))+h_x_pow2*(T_new(i,j-1)+T_new(i,j+1)) ) ...
                )...
            )^2;
    end
end
residual = 1/(N_x*N_y)*sqrt(residual);

end

