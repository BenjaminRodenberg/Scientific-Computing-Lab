close all;
tau=.01;
T_end=10;
t=0:tau:T_end;

tau_ee=1;
tau_h=1;
f0=1;

[t_ee,y_ee]=explicit_euler(@dp,f0,tau_ee,T_end);
[t_h,y_h]=heun(@dp,f0,tau_h,T_end);
[t_ode45,y_ode45]=ode45(@dp,[0 T_end],f0);

figure(1)
hold on
h_ana_p     =plot(t,p(t),'k');
h_ana_dp    =plot(t,dp(t,p(t)),'k:');
h_ee_cont   =plot(t_ee,y_ee,'r');
h_ee_pt     =plot(t_ee,y_ee,'rx');
h_h_cont    =plot(t_h,y_h,'b');
h_h_pt      =plot(t_h,y_h,'bx');
legend([h_ana_p,h_ana_dp,h_ee_pt,h_h_pt],'p(x) analytical','dp(x)','p(x) explicit euler','p(x) heun');
hold off

figure(2)
hold on
title('error')
h_e_ee      =plot(t_ee,abs(y_ee-p(t_ee)),'r');
h_e_h       =plot(t_h,abs(y_h-p(t_h)),'b');
legend([h_e_ee,h_e_h],'error explicit euler','error heun');
hold off