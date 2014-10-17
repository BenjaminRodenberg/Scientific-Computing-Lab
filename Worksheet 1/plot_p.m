close all;
clear all;

euler_index = 1;
heun_index = 2;
runge_kutta_index = 3;

t_end=5;
%t_end=5;

%tau=t_end*1./2.^(10:-1:0); 
tau = [1/8,1/4,1/2,1];
tau_count=numel(tau);

E=struct();

E.compareAna.tau=zeros(tau_count,1);
E.compareAna.euler=zeros(tau_count,1);
E.compareAna.heun=zeros(tau_count,1);
E.compareAna.rungekutta=zeros(tau_count,1);
E.compareNum.tau=zeros(tau_count,1);
E.compareNum.euler=zeros(tau_count,1);
E.compareNum.heun=zeros(tau_count,1);
E.compareNum.rungekutta=zeros(tau_count,1);

n_schemes=3;
y0=1;
y_best = 0;
t_best = 0;

for sim=1:tau_count

curr_tau=tau(sim);

t=0:curr_tau:t_end;

if ~y_best
    y_best = zeros(n_schemes,size(t,2));
    t_best= zeros(n_schemes,size(t,2));
end

[t_ee,y_ee]=explicit_euler(@dp,y0,curr_tau,t_end);
[t_h,y_h]=heun(@dp,y0,curr_tau,t_end);
[t_rk,y_rk]=runge_kutta_4(@dp,y0,curr_tau,t_end);

y_num=zeros(n_schemes,size(t,2));
y_num(euler_index,:)=y_ee;
y_num(heun_index,:)=y_h;
y_num(runge_kutta_index,:)=y_rk;

t_num=zeros(size(y_num));
t_num(euler_index,:)=t_ee;
t_num(heun_index,:)=t_h;
t_num(runge_kutta_index,:)=t_rk;

curr_errors_ana=zeros(n_schemes,1);
for i = 1:n_schemes    
    curr_errors_ana(i)=error_norm(y_num(i,:),p(t_num(i,:)),curr_tau,t_end);
end

E.compareAna.euler(sim)=curr_errors_ana(euler_index);
E.compareAna.heun(sim)=curr_errors_ana(heun_index);
E.compareAna.rungekutta(sim)=curr_errors_ana(runge_kutta_index);

if sim==1 %%Case when curr_tau is minimum
    y_best=y_num;
    t_best=t_num; 
    E.compareNum.euler(sim)=0;         
    E.compareNum.heun(sim)=0;
    E.compareNum.rungekutta(sim)=0;
else   
    y_best_interpolated = zeros( n_schemes, numel(y_num(euler_index,:)));
    y_best_interpolated(euler_index,:) = interp1(t_best(euler_index,:), y_best(euler_index,:), t_num(euler_index,:));
    y_best_interpolated(heun_index,:) = interp1(t_best(heun_index,:), y_best(heun_index,:), t_num(euler_index,:));
    y_best_interpolated(runge_kutta_index,:) = interp1(t_best(runge_kutta_index,:), y_best(runge_kutta_index,:), t_num(euler_index,:));    
    y_best_interpolated(euler_index,:)

    E.compareNum.euler(sim)=error_norm(y_ee,y_best_interpolated(euler_index,:),curr_tau,t_end);         
    E.compareNum.heun(sim)=error_norm(y_h,y_best_interpolated(heun_index,:),curr_tau,t_end);
    E.compareNum.rungekutta(sim)=error_norm(y_rk,y_best_interpolated(runge_kutta_index,:),curr_tau,t_end);
end

end

figure(1)
hold on

title('Error Plot')
xlabel('tau')
ylabel('Error')

h_ee_Eana=plot(tau,E.compareAna.euler,'rx');
plot(tau,E.compareAna.euler,'r');
h_he_Eana=plot(tau,E.compareAna.heun,'bx');
plot(tau,E.compareAna.heun,'b');
h_rk_Eana=plot(tau,E.compareAna.rungekutta,'gx');
plot(tau,E.compareAna.rungekutta,'g');

h_ee_Enum=plot(tau,E.compareNum.euler,'r*');
plot(tau,E.compareNum.euler,'r');
h_he_Enum=plot(tau,E.compareNum.heun,'b*');
plot(tau,E.compareNum.heun,'b');
h_rk_Enum=plot(tau,E.compareNum.rungekutta,'g*');
plot(tau,E.compareNum.rungekutta,'g');

legend([h_ee_Eana,h_he_Eana,h_rk_Eana,h_ee_Enum,h_he_Enum,h_rk_Enum],...
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