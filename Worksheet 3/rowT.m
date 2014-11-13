function [ rowT ] = rowT( i, j, Nx, Ny )
%T Summary of this function goes here
%   Detailed explanation goes here


%T = sin(pi*x)*sin(pi*y);

%Ts surrounding node will be weighted:
%One row has Nx * Ny elements
%Use a vector of zeros and set the four relevant values
Tvector = zeros(Nx*Ny, 1);
%Tset @(i,j) = Tvector((j-1)*Nx+i);

%boundaries not regarded at the moment, as they are 0 in this case

hx = 1/(Nx+1);
hy = 1/(Ny+1);


%set the 5 non-zero vector elements
%vector is adressed as: ((j-1)**Nx+i)

% Tset(i-1,j) = 1/hx^2;
% Tset(i,j) = -2*(hx^2+hy^2)/(hx^2*hy^2);
% Tset(i+1,j) = 1/hx^2;
% Tset(i,j-1) = 1/hy^2;
% Tset(i,j+1) = 1/hy^2;

i_temp=i-1; j_temp=j;
if (i_temp>0 & i_temp<Nx+1 & j_temp>0& j_temp<Ny+1)
    Tvector((j_temp-1)*Nx+i_temp) = 1/hx^2;
end
i_temp=i; j_temp=j;
if (i_temp>0 & i_temp<Nx+1 & j_temp>0& j_temp<Ny+1)
    Tvector((j_temp-1)*Nx+i_temp) = -2*(hx^2+hy^2)/(hx^2*hy^2);
end
i_temp=i+1; j_temp=j;
if (i_temp>0 & i_temp<Nx+1 & j_temp>0& j_temp<Ny+1)
    Tvector((j_temp-1)*Nx+i_temp) = 1/hx^2;
end
i_temp=i; j_temp=j-1;
if (i_temp>0 & i_temp<Nx+1 & j_temp>0& j_temp<Ny+1)
    Tvector((j_temp-1)*Nx+i_temp) = 1/hy^2;
end
i_temp=i; j_temp=j+1;
if (i_temp>0 & i_temp<Nx+1 & j_temp>0& j_temp<Ny+1)
    Tvector((j_temp-1)*Nx+i_temp) = 1/hy^2;
end
% for j =1:1:Ny
%     for i=1:1:Nx
%         Tvector((j-1)*Nx+i) = j;
%     end
% end



Tvector.'
rowT = Tvector;


end

