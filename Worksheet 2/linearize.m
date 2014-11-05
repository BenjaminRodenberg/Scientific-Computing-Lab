function [linearized_f,mat_linearized_f] = linearize(f)

syms x0;
x = symvar( f );
f(x) = f;
Df(x) = gradient( f );
linearized_f( x, x0 ) = simplify( Df(x0)*(x-x0) + f(x0) );
mat_linearized_f = matlabFunction( linearized_f );

end