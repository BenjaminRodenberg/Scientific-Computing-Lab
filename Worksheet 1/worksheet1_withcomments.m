close all;
clear all;

euler_index = 1;
heun_index = 2;
runge_kutta_index = 3;

t_end = 5;

%Different step sizes tau in time domain
tau = [1/8,1/4,1/2,1];
tau_count = numel( tau );

%Structure array to store error values of the different evaluations
E = struct();
%Structure array to store results of the numerical solutions
num_sol = struct();

%Initialize arrays: Error in comparison to exact solution
E.error.euler = zeros( tau_count, 1 );
E.error.heun = zeros( tau_count, 1 );
E.error.rungekutta = zeros( tau_count, 1 );

%Initialize arrays: Approximated error in comparison to best numerical
%solution, using the given error approximation formula
E.approx_error.euler = zeros( tau_count, 1 );
E.approx_error.heun = zeros( tau_count, 1 );
E.approx_error.rungekutta = zeros( tau_count, 1 );

%Initialize arrays: Error reduction for the 3 methods for halfed step size
E.red_error.euler = 0;
E.red_error.heun = 0;
E.red_error.rungekutta = 0;

%Initialize arrays: Error reduction factor with halfed step size in
%comparison to normal step size
E.red_factor.euler = zeros( tau_count, 1 );
E.red_factor.heun = zeros( tau_count, 1 );
E.red_factor.rungekutta = zeros( tau_count, 1 );

%Initialization of number of methods, start value of function and variables
%for the best solution arrays
schemes_count = 3;
y0=1;
y_best = 0;
t_best = 0;
curr_errors_ana = zeros( schemes_count, 1 ); %3x1 array

%Calculation of the results for the different given step sizes
for current_tau_index = 1:tau_count
    curr_tau=tau(current_tau_index);
    %t is the vector which contains all time steps
    t=0:curr_tau:t_end;
    
    %y_best/t_best are arrays with an individual row for each scheme. Each
    %row has as many elements as there are time steps
    if ~y_best
        y_best = zeros(schemes_count,size(t,2));
        t_best= zeros(schemes_count,size(t,2));
    end
    
    %Calculate the approximate values for function dp at the current step
    %size for the three numerical methods and save them in the arrays y_ee,
    %y_h and y_rk for the respective times.
    [t_ee,y_ee] = explicit_euler( @dp, y0, curr_tau, t_end );
    [t_h,y_h] = heun( @dp, y0, curr_tau, t_end );
    [t_rk,y_rk] = runge_kutta_4( @dp, y0, curr_tau, t_end );
    
    %y_num is an array containing the values of the arrays y_ee, y_h and
    %y_rk in its rows.
    y_num = zeros( schemes_count, size(t,2) );
    y_num( euler_index, : ) = y_ee;
    y_num( heun_index, : ) = y_h;
    y_num( runge_kutta_index , : ) = y_rk;
    
    %t_num contains the corresponding time steps to the values in y_num
    t_num = zeros( size( y_num ) );
    t_num( euler_index, : ) = t_ee;
    t_num( heun_index, : ) = t_h;
    t_num( runge_kutta_index, : ) = t_rk;
    
    %Form a string of format tau1, tau2 or tau3 (For tau indices 1 to 3),
    %depending on current tau.
    str_tau=sprintf('tau%i',current_tau_index);
    
    %Save the results of the numerical solution to the structure array
    %num_sol and address them with the current step size (tau).
    num_sol.euler.(str_tau).tau=curr_tau;
    num_sol.euler.(str_tau).y=y_ee;
    num_sol.euler.(str_tau).t=t_ee;
    num_sol.heun.(str_tau).tau=curr_tau;
    num_sol.heun.(str_tau).y=y_h;
    num_sol.heun.(str_tau).t=t_h;
    num_sol.rungekutta.(str_tau).tau=curr_tau;
    num_sol.rungekutta.(str_tau).y=y_rk;
    num_sol.rungekutta.(str_tau).t=t_rk;
    
    %Calculate the errors in comparison to the exact solution for all
    %numerical methods. y_num are the numerically computed function values
    %and p(t) is the exact analytical solution.
    for i = 1:schemes_count
        curr_errors_ana(i) = Error_norm( y_num(i,:), p(t_num(i,:)), curr_tau, t_end);
    end
    
    %Save the calculated errors to the structure array
    E.error.euler( current_tau_index ) = curr_errors_ana( euler_index );
    E.error.heun( current_tau_index ) = curr_errors_ana( heun_index );
    E.error.rungekutta( current_tau_index ) = curr_errors_ana( runge_kutta_index );
    
    %Calculate time and result vectors for the halfed step size
    [t_ee_halved,y_ee_halved] = explicit_euler( @dp, y0, curr_tau/2, t_end );
    [t_h_halved,y_h_halved] = heun( @dp, y0, curr_tau/2, t_end );
    [t_rk_halved,y_rk_halved] = runge_kutta_4( @dp, y0, curr_tau/2, t_end );
    
    %The first if-condition has to be additionally calculated. For the
    %following step sizes the previous results can be used, as tau is
    %doubled with each iteration.
    if current_tau_index == 1 %%Case when curr_tau is minimum
        E.red_error.euler = Error_norm( y_ee_halved , p(t_ee_halved), curr_tau/2, t_end);
        E.red_error.heun = Error_norm(y_h_halved, p(t_h_halved), curr_tau/2, t_end);
        E.red_error.rungekutta = Error_norm( y_rk_halved, p(t_rk_halved), curr_tau/2, t_end);
        
        %Calculate factor between regular error and error with halfed step
        %size. Save results to structure array.
        E.red_factor.euler( current_tau_index ) = E.error.euler( current_tau_index ) / E.red_error.euler;
        E.red_factor.heun( current_tau_index ) = E.error.heun( current_tau_index ) / E.red_error.heun;
        E.red_factor.rungekutta( current_tau_index ) = E.error.rungekutta( current_tau_index ) / E.red_error.rungekutta;
    else
        E.red_factor.euler( current_tau_index ) = E.error.euler( current_tau_index ) / E.error.euler( current_tau_index - 1 );
        E.red_factor.heun( current_tau_index ) = E.error.heun( current_tau_index ) / E.error.heun( current_tau_index - 1 );
        E.red_factor.rungekutta( current_tau_index ) = E.error.rungekutta( current_tau_index ) / E.error.rungekutta( current_tau_index - 1 );
    end
    
    %Use the function/time values of numerical solution with minimum step
    %size as the best approximation.
    if current_tau_index == 1 %%Case when curr_tau is minimum
        y_best=y_num;
        t_best=t_num;
        %No improvement to best numerical solution is possible, as this is
        %the best solution
        E.approx_error.euler(current_tau_index)=0;
        E.approx_error.heun(current_tau_index)=0;
        E.approx_error.rungekutta(current_tau_index)=0;
    else
        %Array that contains as many entries as there are calculated
        %function values for each of the schemes.
        y_best_interpolated = zeros( schemes_count, numel( y_num( euler_index, : ) ) );
        
        %y_best was calculated for many time steps, using the lowest step
        %size. In order to obtain the most appropriate value of the y_best
        %graph at time t_num, y_best is interpolated to provide that value.
        for i=1:schemes_count
            y_best_interpolated(i,:) = interp1( t_best(i,:), y_best(i,:), t_num(i,:));
        end
        
        %Save error in comparison to best numerical solution in structure
        %array
        E.approx_error.euler(current_tau_index) = Error_norm( y_ee,y_best_interpolated(euler_index,:),...
            curr_tau,t_end);
        E.approx_error.heun(current_tau_index) = Error_norm(y_h,y_best_interpolated(heun_index,:),...
            curr_tau,t_end);
        E.approx_error.rungekutta(current_tau_index) = Error_norm(y_rk,y_best_interpolated(runge_kutta_index,:),...
            curr_tau,t_end);
    end
    
end

%Printing of the results

fprintf('Error comparing to exact solution:\n')
fprintf('----------------------------------------------------------------------------\n')
fprintf('tau              |%d\t\t\t\t%d\t%d\t%d\n',tau(end:-1:1))
fprintf('explicit euler   |%d\t%d\t%d\t%d\n',E.error.euler(end:-1:1))
fprintf('heun             |%d\t%d\t%d\t%d\n',E.error.heun(end:-1:1))
fprintf('rungekutta       |%d\t%d\t%d\t%d\n',E.error.rungekutta(end:-1:1))
fprintf('----------------------------------------------------------------------------\n\n')

fprintf('Error reduction:\n')
fprintf('----------------------------------------------------------------------------\n')
fprintf('tau              |%d\t\t\t\t%d\t%d\t%d\n',tau(end:-1:1))
fprintf('explicit euler   |%d\t%d\t%d\t%d\n',E.red_factor.euler(end:-1:1))
fprintf('heun             |%d\t%d\t%d\t%d\n',E.red_factor.heun(end:-1:1))
fprintf('rungekutta       |%d\t%d\t%d\t%d\n',E.red_factor.rungekutta(end:-1:1))
fprintf('----------------------------------------------------------------------------\n\n')

fprintf('Error comparing to best numerical solution:\n')
fprintf('----------------------------------------------------------------------------\n')
fprintf('tau              |%d\t\t\t\t%d\t%d\t%d\n',tau(end:-1:1))
fprintf('explicit euler   |%d\t%d\t%d\t%d\n',E.approx_error.euler(end:-1:1))
fprintf('heun             |%d\t%d\t%d\t%d\n',E.approx_error.heun(end:-1:1))
fprintf('rungekutta       |%d\t%d\t%d\t%d\n',E.approx_error.rungekutta(end:-1:1))
fprintf('----------------------------------------------------------------------------\n\n')

%plot of analytical p(t)
figure(1)
hold on
title('p(t)')
plot(0:.1:t_end,p(0:.1:t_end),'k')
xlabel('t')
ylabel('p(t)')
hold off

%plot of numerical solutions using explicit euler method
figure(2)
hold on
title('numerical solution of dp=(1-p/10)*p using explicit euler method');
h_ana=plot(0:.1:t_end,p(0:.1:t_end),'k');
col={'r.','rx','r*','rs'};
for i = 1:4    
    field_name=sprintf('tau%i',i);
    %The values for the numerical calculation with different step sizes
    %have been saved in the num_sol structure array. They are adressed
    %individually in each iteration by using the strings tau1 to tau4.
    h(i)=plot(num_sol.euler.(field_name).t,num_sol.euler.(field_name).y,col{i});
end
legend([h_ana,h],'analytical solution',...
    ['numerical solution for tau=',mat2str(num_sol.euler.tau1.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.euler.tau2.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.euler.tau3.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.euler.tau4.tau)]);
xlabel('t')
ylabel('p(t)')
hold off

%plot of numerical solutions using heun method
figure(3)
hold on
title('numerical solution of dp=(1-p/10)*p using heun method');
h_ana=plot(0:.1:t_end,p(0:.1:t_end),'k');
col={'b.','bx','b*','bs'};
for i = 1:4    
    field_name=sprintf('tau%i',i);
    h(i)=plot(num_sol.heun.(field_name).t,num_sol.heun.(field_name).y,col{i});
end
legend([h_ana,h],'analytical solution',...
    ['numerical solution for tau=',mat2str(num_sol.heun.tau1.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.heun.tau2.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.heun.tau3.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.heun.tau4.tau)]);
xlabel('t')
ylabel('p(t)')
hold off

%plot of numerical solutions using runge-kutta method
figure(4)
hold on
title('numerical solution of dp=(1-p/10)*p using runge-kutta method');
h_ana=plot(0:.1:t_end,p(0:.1:t_end),'k');
col={'g.','gx','g*','gs'};
for i = 1:4    
    field_name=sprintf('tau%i',i);
    h(i)=plot(num_sol.rungekutta.(field_name).t,num_sol.rungekutta.(field_name).y,col{i});
end
legend([h_ana,h],'analytical solution',...
    ['numerical solution for tau=',mat2str(num_sol.rungekutta.tau1.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.rungekutta.tau2.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.rungekutta.tau3.tau)],...
    ['numerical solution for tau=',mat2str(num_sol.rungekutta.tau4.tau)]);
xlabel('t')
ylabel('p(t)')
hold off

%comparison of the error
figure(5)
hold on

title('Error Plot')
xlabel('tau')
ylabel('Error')

h_euler_error = plot( tau, E.error.euler, 'rx' );
plot( tau, E.error.euler, 'r' );
h_heun_error = plot( tau,E.error.heun,'bx' );
plot( tau, E.error.heun, 'b' );
h_rungekutta_error = plot( tau, E.error.rungekutta, 'gx' );
plot( tau, E.error.rungekutta, 'g');

h_euler_approx_error = plot( tau, E.approx_error.euler, 'r*' );
plot( tau, E.approx_error.euler, 'r' );
h_heun_approx_error = plot( tau, E.approx_error.heun, 'b*' );
plot( tau,E.approx_error.heun,'b');
h_rungekutta_approx_error=plot(tau,E.approx_error.rungekutta,'g*');
plot( tau,E.approx_error.rungekutta,'g');

legend([h_euler_error,h_heun_error,h_rungekutta_error,h_euler_approx_error,...
        h_heun_approx_error,h_rungekutta_approx_error],...
        'Error explicit euler analytical',...
        'Error heun analytical',...
        'Error runge kutta analytical',...
        'Error explicit euler numerical',...
        'Error heun numerical',...
        'Error runge kutta numerical');

% set(gca,'Xscale','log')
% set(gca,'Yscale','log')

hold off