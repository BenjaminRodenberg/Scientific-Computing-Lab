function [ x ] = newton( f,df,x0,tol )
%NEWTON newton method for solving root equation numerically

x=x0;
%error at the beginning = inf
errornew=inf;

while(errornew>tol)
    
    xold=x;
    %newton method fixpoint equation
    x=x-f(x)/df(x);
    
    %calculate error via difference of approximate solution
    errorold=errornew;
    errornew=abs(x-xold);
    
    %check weather newton method does converge. If newton method does
    %converge, the error decreases monotonically.
    if(errornew>errorold)
        errornew('newton iteration does not converge!');
    end
        
end

end

