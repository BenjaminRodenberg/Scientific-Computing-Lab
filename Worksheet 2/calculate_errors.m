function [ sols ] = calculate_errors( sols )
%CALCULATE_ERRORS calculates the L2-Error of numerical solutions of an ODE
%using different ODE solvers. Finally also the error reduction with
%decreasing timestep size is calculated.

errors = calculate_analitical_errors( sols );

%save results
sols.errors = errors;

%save result
sols.error_reduction = calculate_error_reduction( errors, sols.tau );

end
