function [ T, residual ] = do_one_Gauss_Seidl_Iteration( N_x,N_y,b,T,BC, h_x,h_y )
%DO_GAUSS_SEIDL_ITERATION Calculates the result of one Gauss Seidl
%Iteration. To save computation time the Residual is calculated on the fly.
%But with a delay of two rows so that neighbour values do not change
%anymore during the current iteration step.

N=(N_x+2*BC)*(N_y+2*BC);

residual = 0;
for i = 1:N_x + 2*BC
    for j = 1:N_y + 2*BC
        %get index in linear array from two-D coordinates i,j
        current_node_index = i + (N_x + 2*BC) * (j-1);
        
        %local residual set to zero: 0 = b-weights*neighbouring Temperature
        %This can be solved for the temperature on the current node.
        
        %S is the sum of the weights times the Temperatures on the
        %neighbour nodes.
        S=0;

        %decide weather node is next to boundary node
        x_up = i<N_x + 2*BC;
        x_dwn = i>1;
        y_up = j<N_y + 2*BC;
        y_dwn = j>1;
        
        %add summands due to boundary conditions
        if(x_up)
            S=S+1/h_x^2*T(current_node_index+1);
        end
        if(x_dwn)
            S=S+1/h_x^2*T(current_node_index-1);
        end
        if(y_up)
            S=S+1/h_y^2*T(current_node_index+N_x);
        end
        if(y_dwn)
            S=S+1/h_y^2*T(current_node_index-N_x);
        end
        
        %solve for T_ij
        T(current_node_index) = (b(current_node_index)-S) * (h_x^2*h_y^2/(-2*(h_x^2+h_y^2)));
        
        %Residual: go two rows back
        if j>2
            %delay of two rows
            j_displaced = j-2;
            current_node_index = i + (N_x + 2*BC) * (j_displaced-1);
            
            S=0;

            %decide weather node is next to boundary node
            x_up = i<N_x + 2*BC;
            x_dwn = i>1;
            y_up = j_displaced<N_y + 2*BC;
            y_dwn = j_displaced>1;
            
            %add summands due to boundary conditions
            if(x_up)
                S=S+1/h_x^2*T(current_node_index+1);
            end
            if(x_dwn)
                S=S+1/h_x^2*T(current_node_index-1);
            end
            if(y_up)
                S=S+1/h_y^2*T(current_node_index+N_x);
            end
            if(y_dwn)
                S=S+1/h_y^2*T(current_node_index-N_x);
            end
            %add weight*T_ij to the sum
            S=S-2*(h_x^2+h_y^2)/(h_x^2*h_y^2)*T(current_node_index);
            
            %iteratively sum up residual
            residual = residual + (b(current_node_index)-S)^2;
            
        end
    end
end

%calculate remaining two rows
for i = 1:N_x + 2*BC
    for j = N_y-1:N_y + 2*BC
        current_node_index = i + (N_x + 2*BC) * (j-1);
        
        S=0;
        
        %decide weather node is next to boundary node
        x_up = i<N_x + 2*BC;
        x_dwn = i>1;
        y_up = j<N_y + 2*BC;
        y_dwn = j>1;
        
        %add summands due to boundary conditions
        if(x_up)
            S=S+1/h_x^2*T(current_node_index+1);
        end
        if(x_dwn)
            S=S+1/h_x^2*T(current_node_index-1);
        end
        if(y_up)
            S=S+1/h_y^2*T(current_node_index+N_x);
        end
        if(y_dwn)
            S=S+1/h_y^2*T(current_node_index-N_x);
        end
        
        S=S-2*(h_x^2+h_y^2)/(h_x^2*h_y^2)*T(current_node_index);
        
        residual = residual + (b(current_node_index)-S)^2;
    end
end

%normalize residual
residual=sqrt(1/N*residual);

return

function [S]=get_sum_neighbours(T,N_x,N_y,h_x,h_y,i,j,current_node_index,BC)

S=0;

%decide weather node is next to boundary node
x_up = i<N_x + 2*BC;
x_dwn = i>1;
y_up = j<N_y + 2*BC;
y_dwn = j>1;

%add summands due to boundary conditions
if(x_up)
    S=S+1/h_x^2*T(current_node_index+1);
end
if(x_dwn)
    S=S+1/h_x^2*T(current_node_index-1);
end
if(y_up)
    S=S+1/h_y^2*T(current_node_index+N_x);
end
if(y_dwn)
    S=S+1/h_y^2*T(current_node_index-N_x);
end

return

