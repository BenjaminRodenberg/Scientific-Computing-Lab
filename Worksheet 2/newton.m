function [ new_x,x_history ] = newton( f,df,x0,tol,calling_function )
%NEWTON newton method for solving root equation numerically

new_x = x0;
%error at the beginning = inf
curr_error = inf;
x_history = new_x;

%newton method fixpoint equation
phi=@(x)x-f(x)/df(x);

while( curr_error > tol && df(new_x)~=0)
    
    old_x = new_x;      
    
    new_x = phi(old_x);    
        
    x_history = [x_history new_x];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The usual convergence test has the following form : 
    %
    % if( abs(new_f_x)>abs(old_f_x) ) STOP
    %
    % This criterion does not lead to a correct result for the following
    % input, even if the method converges:
    %
    % newton(@(x)1-x^2,@(x)-2*x,.1,10^-4,'')
    % 
    % Therefore we are using the contraction condition from banach's
    % fixpoint theorem. The following examples show the correct behaviour:
    %    
    % newton(@(x)atan(x),@(x)1/(x^2+1),.4,10^-4,'') 
    % newton(@(x)atan(x),@(x)1/(x^2+1),1.6,10^-4,'') 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % check weather newton method does converge. 
    if( abs(phi(new_x)-phi(old_x))/abs(new_x-old_x) > 1 )
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

