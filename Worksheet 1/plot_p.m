close all;
clear all;

euler_index = 1;
heun_index = 2;
runge_kutta_index = 3;

t_end = 5;

tau = [1/8,1/4,1/2,1];
tau_count = numel( tau );

E = struct();

E.error.euler = zeros( tau_count, 1 );
E.error.heun = zeros( tau_count, 1 );
E.error.rungekutta = zeros( tau_count, 1 );

E.approx_error.euler = zeros( tau_count, 1 );
E.approx_error.heun = zeros( tau_count, 1 );
E.approx_error.rungekutta = zeros( tau_count, 1 );

schemes_count = 3;
y0=1;
y_best = 0;
t_best = 0;
curr_errors_ana = zeros( schemes_count, 1 );

for current_tau_index = 1:tau_count
curr_tau=tau(current_tau_index);
t=0:curr_tau:t_end;

if ~y_best
    y_best = zeros(schemes_count,size(t,2));
    t_best= zeros(schemes_count,size(t,2));
end

[t_ee,y_ee] = explicit_euler( @dp, y0, curr_tau, t_end );
[t_h,y_h] = heun( @dp, y0, curr_tau, t_end );
[t_rk,y_rk] = runge_kutta_4( @dp, y0, curr_tau, t_end );

y_num = zeros( schemes_count, size(t,2) );
y_num( euler_index, : ) = y_ee;
y_num( heun_index, : ) = y_h;
y_num( runge_kutta_index , : ) = y_rk;

t_num = zeros( size( y_num ) );
t_num( euler_index, : ) = t_ee;
t_num( heun_index, : ) = t_h;
t_num( runge_kutta_index, : ) = t_rk;

for i = 1:schemes_count    
    curr_errors_ana(i) = error_norm( y_num(i,:), p(t_num(i,:)), curr_tau, t_end);
end

E.error.euler( current_tau_index ) = curr_errors_ana( euler_index );
E.error.heun( current_tau_index ) = curr_errors_ana( heun_index );
E.error.rungekutta( current_tau_index ) = curr_errors_ana( runge_kutta_index );

[t_ee_halved,y_ee_halved] = explicit_euler( @dp, y0, curr_tau/2, t_end );
[t_h_halved,y_h_halved] = heun( @dp, y0, curr_tau/2, t_end );
[t_rk_halved,y_rk_halved] = runge_kutta_4( @dp, y0, curr_tau/2, t_end );

E.red_error.euler( current_tau_index ) = error_norm( y_ee_halved(1:2:numel(y_ee_halved)), p(t_num(euler_index,:)), curr_tau, t_end);
E.red_error.heun( current_tau_index ) = error_norm(y_h_halved(1:2:numel(y_ee_halved)), p(t_num(heun_index,:)), curr_tau, t_end);
E.red_error.rungekutta( current_tau_index ) = error_norm( y_rk_halved(1:2:numel(y_ee_halved)), p(t_num(runge_kutta_index,:)), curr_tau, t_end);

if current_tau_index == 1 %%Case when curr_tau is minimum
    y_best=y_num;
    t_best=t_num; 
    E.approx_error.euler(current_tau_index)=0;         
    E.approx_error.heun(current_tau_index)=0;
    E.approx_error.rungekutta(current_tau_index)=0;
else   
    y_best_interpolated = zeros( schemes_count, numel( y_num( euler_index, : ) ) );
    
    for i=1:schemes_count
        y_best_interpolated(i,:) = interp1( t_best(i,:), y_best(i,:), t_num(i,:));
    end

    E.approx_error.euler(current_tau_index) = error_norm( y_ee,y_best_interpolated(euler_index,:),...
                                                        curr_tau,t_end);         
    E.approx_error.heun(current_tau_index) = error_norm(y_h,y_best_interpolated(heun_index,:),...
                                                        curr_tau,t_end);
    E.approx_error.rungekutta(current_tau_index) = error_norm(y_rk,y_best_interpolated(runge_kutta_index,:),...
                                                        curr_tau,t_end);
end

end

figure(1)
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

h_euler_red_error = plot( tau, E.red_error.euler, 'ro-' );
plot( tau, E.red_error.euler, 'r' );
h_heun_red_error = plot( tau, E.red_error.heun, 'bo-' );
plot( tau,E.red_error.heun,'b');
h_rungekutta_red_error=plot(tau,E.red_error.rungekutta,'go-');
plot( tau,E.red_error.rungekutta,'g');

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

% figure(1)
% hold on
% h_ana_p     =plot(t,p(t),'k');
% h_ana_dp    =plot(t,dp(t,p(t)),'k:');
% h_ee_cont   =plot(t_ee,y_ee,'r');
% h_ee_pt     =plot(t_ee,y_ee,'rx');
% h_h_cont    =plot(t_h,y_h,'b');
% h_h_pt      =plot(t_h,y_h,'bx');
% h_rk_cont   =plot(t_rk,y_rk,'g');
% h_rk_pt     =plot(t_rk,y_rk,'gx');
% legend([h_ana_p,h_ana_dp,h_ee_pt,h_h_pt,h_rk_pt],'p(x) analytical','dp(x)','p(x) explicit euler','p(x) heun','p(x) runge kutta 4th');
% hold off
% 
% figure(2)
% hold on
% title('error')
% h_e_ee      =plot(t_ee,abs(y_ee-p(t_ee)),'r');
% h_e_h       =plot(t_h,abs(y_h-p(t_h)),'b');
% h_e_rk      =plot(t_h,abs(y_rk-p(t_rk)),'g');
% legend([h_e_ee,h_e_h,h_e_rk],'error explicit euler','error heun','error runge kutta 4th');
% hold off