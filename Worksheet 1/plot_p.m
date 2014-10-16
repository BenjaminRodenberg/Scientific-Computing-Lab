close all;
%T_end=10;
T_end=5;
curr_tau=T_end/100;
t=0:curr_tau:T_end;

%tau_num=T_end*1./2.^(10:-1:0);
tau_num=[1/8,1/4,1/2,1];
n_sim=numel(tau_num);

E=struct();
E.compareAna.tau=zeros(n_sim,1);
E.compareAna.euler=zeros(n_sim,1);
E.compareAna.heun=zeros(n_sim,1);
E.compareAna.rungekutta=zeros(n_sim,1);
E.compareNum.tau=zeros(n_sim,1);
E.compareNum.euler=zeros(n_sim,1);
E.compareNum.heun=zeros(n_sim,1);
E.compareNum.rungekutta=zeros(n_sim,1);

for sim=1:n_sim
    
    curr_tau=tau_num(sim);
    
    y0=1;
    
    n_schemes=3;
    
    t=0:curr_tau:T_end;
    
    [t_ee,y_ee]=explicit_euler(@dp,y0,curr_tau,T_end);
    [t_h,y_h]=heun(@dp,y0,curr_tau,T_end);
    [t_rk,y_rk]=runge_kutta_4(@dp,y0,curr_tau,T_end);
    
    %[t_ode45,y_ode45]=ode45(@dp,[0 T_end],y0);
    
    y_num=zeros(n_schemes,size(t,2));
    y_num(1,:)=y_ee;
    y_num(2,:)=y_h;
    y_num(3,:)=y_rk;
    
    t_num=zeros(n_schemes,size(t,2));
    t_num(1,:)=t_ee;
    t_num(2,:)=t_h;
    t_num(3,:)=t_rk;
    
    curr_Errors_ana=zeros(n_schemes,1);
    for i = 1:n_schemes
        curr_Errors_ana(i)=Error_norm(y_num(i,:),p(t_num(i,:)),curr_tau,T_end);
    end
    
    E.compareAna.tau(sim)=curr_tau;
    E.compareAna.euler(sim)=curr_Errors_ana(1);
    E.compareAna.heun(sim)=curr_Errors_ana(2);
    E.compareAna.rungekutta(sim)=curr_Errors_ana(3);

    fieldname=['tau',mat2str(sim)];
    sol.euler.(fieldname).y=y_num(1,:);
    sol.euler.(fieldname).t=t_num(1,:);
    sol.heun.(fieldname).y=y_num(2,:);
    sol.heun.(fieldname).t=t_num(2,:);
    sol.rungekutta.(fieldname).y=y_num(3,:);
    sol.rungekutta.(fieldname).t=t_num(3,:);
        
end

y_best=zeros(n_schemes,floor(T_end/min(tau_num))+1);

best_sim=find(min(tau_num)==tau_num);
fieldname=['tau',num2str(best_sim)];    
y_best(1,:)=sol.euler.(fieldname).y;
t_best(1,:)=sol.euler.(fieldname).t;
y_best(2,:)=sol.heun.(fieldname).y;
t_best(2,:)=sol.heun.(fieldname).t;
y_best(3,:)=sol.rungekutta.(fieldname).y;
t_best(3,:)=sol.rungekutta.(fieldname).t;

for sim = 1:n_sim  
     
    curr_tau=tau_num(sim);

    fieldname=['tau',num2str(sim)];    

    y_ee=sol.euler.(fieldname).y;
    t_ee=sol.euler.(fieldname).t;
    y_he=sol.heun.(fieldname).y;
    t_he=sol.heun.(fieldname).t;
    y_rk=sol.rungekutta.(fieldname).y;
    t_rk=sol.rungekutta.(fieldname).t;
    
    t_best_ee=t_best(1,:);
    y_best_ee=y_best(1,:);
    y_best_ee=y_best_ee(arrayfun(@(x) find(t_best_ee == x,1,'first'), t_ee ));
    t_best_he=t_best(2,:);
    y_best_he=y_best(2,:);
    y_best_he=y_best_he(arrayfun(@(x) find(t_best_he == x,1,'first'), t_he ));
    t_best_rk=t_best(3,:);
    y_best_rk=y_best(3,:);
    y_best_rk=y_best_rk(arrayfun(@(x) find(t_best_rk == x,1,'first'), t_rk ));

    E.compareNum.tau(sim)       =curr_tau;
    E.compareNum.euler(sim)     =Error_norm(y_ee,y_best_ee,curr_tau,T_end);
    E.compareNum.heun(sim)      =Error_norm(y_he,y_best_he,curr_tau,T_end);
    E.compareNum.rungekutta(sim)=Error_norm(y_rk,y_best_rk,curr_tau,T_end);    
    
end

figure(1)
hold on

title('Error Plot')
xlabel('tau')
ylabel('Error')

h_ee_Eana=plot(E.compareAna.tau,E.compareAna.euler,'rx');
plot(E.compareAna.tau,E.compareAna.euler,'r');
h_he_Eana=plot(E.compareAna.tau,E.compareAna.heun,'bx');
plot(E.compareAna.tau,E.compareAna.heun,'b');
h_rk_Eana=plot(E.compareAna.tau,E.compareAna.rungekutta,'gx');
plot(E.compareAna.tau,E.compareAna.rungekutta,'g');

h_ee_Enum=plot(E.compareNum.tau,E.compareNum.euler,'r*');
plot(E.compareNum.tau,E.compareNum.euler,'r');
h_he_Enum=plot(E.compareNum.tau,E.compareNum.heun,'b*');
plot(E.compareNum.tau,E.compareNum.heun,'b');
h_rk_Enum=plot(E.compareNum.tau,E.compareNum.rungekutta,'g*');
plot(E.compareNum.tau,E.compareNum.rungekutta,'g');

legend([h_ee_Eana,h_he_Eana,h_rk_Eana,h_ee_Enum,h_he_Enum,h_rk_Enum],...
    'Error explicit euler analytical',...
    'Error heun analytical',...
    'Error runge kutta analytical',...
    'Error explicit euler numerical',...
    'Error heun numerical',...
    'Error runge kutta numerical');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')

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