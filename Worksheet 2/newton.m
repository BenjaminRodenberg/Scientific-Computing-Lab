function [ x, x_history ] = newton( f,df,x0,tol )
%NEWTON newton method for solving root equation numerically

x=x0;
%error at the beginning = inf
new_error = inf;
x_history = x;

while( new_error > tol )
    
    old_x = x;
    %newton method fixpoint equation
    x = x - ( f(x) / df(x) );
    
    %calculate error via difference of approximate solution
    old_error = new_error;
    new_error = abs( x - old_x );
    x_history = [x_history x];
    
    %check weather newton method does converge. If newton method does
    %converge, the error decreases monotonically.
    if( new_error >= old_error )
        warning('newton iteration does not converge!');
        break;
    end
        
end

end

