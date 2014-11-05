function [ G,dG ] = derive_G_implicit_euler( sym_f, tau )
%calculates the necessary functions G and dG for reformulating an ODE as a
%fix-point iteration using the newton method.

syms yn yn1
x = symvar( sym_f ); %arguments of RHS
sym_f(x)=sym_f;

sym_G( yn1, yn) = yn + tau * sym_f(yn1) - yn1;
sym_dG( yn1, yn) = diff( sym_G, yn1 );

G = matlabFunction( simplify( sym_G ) );
dG = matlabFunction( simplify( sym_dG ) );

end

