function [ G,dG ] = derive_G_implicit_euler( sym_f,tau )
%calculates the nevessary functions G and dG for reformulating an ODE as a
%fix-point iteration using the newton method.

syms yn
yn1=symvar(sym_f); %arguments of RHS
sym_f(yn1)=sym_f;  %declare symbolic expression f

% root point equation
% G(y_{n+1},y_n) = 0 = y_n+tau*f(y_{n+1})-y_{n+1}
sym_G(yn1,yn)=yn+tau*sym_f(yn1)-yn1;

% differentiate G(y_{n+1},y_n)
sym_dG(yn1,yn)=diff(sym_G,yn1);

% convert symbolic functions into matlab functions
G=matlabFunction(simplify(sym_G));
dG=matlabFunction(simplify(sym_dG));

end

