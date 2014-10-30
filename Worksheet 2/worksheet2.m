%symbolic formulation of the ODE dp(p)=7*(1-p/10)*p
syms sym_p;
sym_dp(sym_p)=7*(1-sym_p/10)*sym_p;

%conversion from symbolic form to matlab function form:
%i.e.: dp=@(p)7*(1-p/10).*p;
dp=matlabFunction(sym_dp);
p0=20;
T_end=5;

%analytical solution of ODE
p_ana=@(t)200./(20-10*exp(-7*t));
t_ana=0:.001:T_end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) Plot the function p(t) in a graph

figure(1)
hold on
title('exact solution of ODE dp(p)=7*(1-p/10)*p')
plot(t_ana,p_ana(t_ana));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

explicit_euler_index    = 1;
heun_index              = 2;
implicit_euler_index    = 3;
adams_moulton_index     = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b) Reuse the Euler method and the method of Heun ...

%different timestep size
tau_range=(.5).^(0:1:5);

numerical_solutions.explicit_euler.tau=tau_range;
for i = 1:numel(tau_range)
    tau=tau_range(i);
    [t_ee,y_ee]=explicit_euler(dp,p0,tau,T_end);
    [t_he,y_he]=heun(dp,p0,tau,T_end);
    field_name=sprintf('tau%i',i);
    numerical_solutions.explicit_euler.(field_name).y=y_ee;
    numerical_solutions.explicit_euler.(field_name).t=t_ee; 
    numerical_solutions.explicit_euler.(field_name).tau=tau;
    numerical_solutions.heun.(field_name).y=y_he;
    numerical_solutions.heun.(field_name).t=t_he; 
    numerical_solutions.heun.(field_name).tau=tau;
end


figure(2)
hold on
title('numerical solution of dp=7*(1-p/10)*p using explicit euler method');

h(1)=plot(t_ana,p_ana(t_ana),'k');
leg_entry{1}='analytical solution';

%col={'r.','rx','r*','rs','r^','rv'};

for i = 1:numel(tau_range)
    field_name=sprintf('tau%i',i);
    h(i+1)=plot(numerical_solutions.explicit_euler.(field_name).t,...
        numerical_solutions.explicit_euler.(field_name).y);%,col{i});
    leg_entry{i+1}=['numerical solution for tau=',...
        mat2str(numerical_solutions.explicit_euler.(field_name).tau)];
    c=[(i-1)/(numel(tau_range)-1) 1-(i-1)/(numel(tau_range)-1) 0];
    set(h(i+1),'Color',c);
end
legend(h,leg_entry);
xlabel('t')
ylabel('p(t)')
axis([0 5 0 20]);
hold off

figure(3)
hold on
title('numerical solution of dp=7*(1-p/10)*p using heun method');

h(1)=plot(t_ana,p_ana(t_ana),'k');
leg_entry{1}='analytical solution';

%col={'b.','bx','b*','bs','b^','bv'};

for i = 1:numel(tau_range)
    field_name=sprintf('tau%i',i);
    h(i+1)=plot(numerical_solutions.heun.(field_name).t,...
        numerical_solutions.heun.(field_name).y);%,col{i});
    leg_entry{i+1}=['numerical solution for tau=',...
        mat2str(numerical_solutions.heun.(field_name).tau)];
    c=[(i-1)/(numel(tau_range)-1) 1-(i-1)/(numel(tau_range)-1) 0];
    set(h(i+1),'Color',c);
end
legend(h,leg_entry);
xlabel('t')
ylabel('p(t)')
axis([0 5 0 20]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [t_exp,y_exp]=explicit_euler(dp,y0,tau,T_end);
% [t_imp,y_imp]=implicit_euler(sym_dp,y0,tau,T_end);
% 
% figure(1)
% subplot(2,1,1)
% hold on
% title('explicit');
% plot(t_exp,y_exp,'b--')
% plot(t_ana,y_ref,'g')
% legend('ODE','reference');
% hold off
% 
% subplot(2,1,2)
% hold on
% title('implicit');
% plot(t_imp,y_imp,'r--')
% plot(t_ana,y_ref,'g')
% legend('ODE','reference');
% hold off
% 
% hold off