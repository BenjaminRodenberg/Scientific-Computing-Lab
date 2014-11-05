function [ new_x,x_history ] = newton( f,df,x0,tol,calling_function )
%NEWTON newton method for solving root equation numerically

new_x = x0;
%error at the beginning = inf
curr_error = inf;
x_history = new_x;

new_f_x = f(new_x);
new_df_x = df(new_x);

while( curr_error > tol && new_df_x~=0)
    
    old_x = new_x;
    old_f_x = new_f_x;
    old_df_x = new_df_x;
    
    %newton method fixpoint equation
    new_x = old_x - (  old_f_x / old_df_x );                    
    new_f_x = f(new_x);
    new_df_x=df(new_x);
    
    x_history = [x_history new_x];
    
    %check weather newton method does converge. If newton method does
    %converge, the error decreases monotonically.
    if( abs(new_f_x) > abs(old_f_x) )
        %calculate error via difference of approximate solution
        curr_error = abs( new_x - old_x );
        warning(['for ',calling_function,' at x = ',mat2str(new_x),...
            ' newton iteration does not converge! Given Tolerance of ',...
            mat2str(tol),' is not met. Error at this point: ',...
            mat2str(curr_error)]);
        
        %return old_x since it is better than new_x
        new_x=old_x;
        break;
    end
    
    %calculate error via difference of approximate solution
    curr_error = abs( new_x - old_x );
            
end

end

