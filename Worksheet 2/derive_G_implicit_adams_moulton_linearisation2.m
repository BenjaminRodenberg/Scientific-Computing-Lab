function [ G,dG ] = derive_G_implicit_adams_moulton_linearisation2( sym_f, tau )
%calculates the necessary functions G and dG for reformulating an ODE as a
%fix-point iteration using the newton method.

syms yn yn1
x = symvar( sym_f ); %arguments of RHS
sym_f(x)=sym_f;

% taking linearisation of adams moulton 
% for f(x)=x*(1-x) the expression below:
% f(yn,yn1)=simplify(sym_f(yn)/yn)*yn1 )
% creates f(yn,yn1)=yn1*(1-yn)
sym_G( yn1, yn ) = yn + (tau/2) * ( sym_f(yn) + simplify(sym_f(yn)/yn)*yn1 ) - yn1;
sym_dG( yn1, yn ) = diff( sym_G, yn1);

G = matlabFunction( simplify( sym_G ) );
dG = matlabFunction( simplify( sym_dG ) );

end
