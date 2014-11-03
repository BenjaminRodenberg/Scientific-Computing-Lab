function [ x,history ] = newton( f,df,x0,tol )
%NEWTON newton method for solving root equation numerically

x=x0;
%error at the beginning = inf
errornew=inf;
history=x;

while(errornew>tol)
    
    xold=x;
    %newton method fixpoint equation
    x=x-f(x)/df(x);
    
    %calculate error via difference of approximate solution
    errorold=errornew;
    errornew=abs(x-xold);
    history=[history x];
    
    %check whether newton method does converge. If newton method does
    %converge, the error decreases monotonically.
    if(errornew>=errorold)
        warning('newton iteration does not converge!');
        break;
    end
        
end

end

