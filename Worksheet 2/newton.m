function [ x ] = newton( f,df,x0,tol )
%NEWTON Summary of this function goes here
%   Detailed explanation goes here

xold=inf;
x=x0;
%history=x;

while(abs(x-xold)>tol)
    
    xold=x;
    x=x-f(x)/df(x);
    %history=[history,x];
    
end

end

