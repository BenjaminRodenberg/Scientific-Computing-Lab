function [ error_reductions ] = calculate_error_reduction( errors, tau_range )

%calculate error reduction
error_reductions = zeros( 1,numel(tau_range) );
error_reductions(1) = nan;

for i = 2:( numel(tau_range) )
    if( errors(i) == inf || errors(i-1) == inf)
        %for inf error no reduction can be calculated
        error_reductions(i) = nan;
    else
        %reduction factor is the ration of old error and new error. If the
        %error is reduced with decreasing timestep size the reduction
        %factor is greater 1.
        error_reductions(i) = errors(i-1)/errors(i);
    end
end