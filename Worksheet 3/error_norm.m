function [E] = error_norm(T_num,T_ref)
%Error norm considering the error of a numerical solution by comparison to
%a better numerical solution or (if available) to the exact analytical
%solution.

N=numel(T_num);

E=sqrt(1/N*sum(sum((T_num-T_ref).^2)));

end