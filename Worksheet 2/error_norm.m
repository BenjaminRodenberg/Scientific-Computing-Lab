function [E] = error_norm(y_num,y_ref,tau,T_end)
%error_norm considering the error of a numerical solution by comparison to
%a better numerical solution or (if available) to the exact analytical
%solution.

E=sqrt(tau/T_end*sum((y_num-y_ref).^2));

end