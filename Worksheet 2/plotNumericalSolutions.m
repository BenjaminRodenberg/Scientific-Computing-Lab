function [ fig_h ] = plotNumericalSolutions( numerical_solutions,method_name,fig_no)
%PLOTNUMERICALSOLUTION plots numerical solutions using the same solver for
%different timestep sizes into figure(fig_no) and compares the numerical
%solution to analytical solution

%analytical solution
y_ana=numerical_solutions.(method_name).y_ana;
t_ana=numerical_solutions.(method_name).t_ana;

%considered timestep sizes
tau_range=numerical_solutions.(method_name).tau;

%open figure
fig_h=figure(fig_no);
hold on

title(['numerical solution of dp=7*(1-p/10)*p using ',method_name,' method'],...
    'Interpreter','none');

%plot analytical solution
h(1)=plot(t_ana,y_ana,'k');
leg_entry{1}='analytical solution';

%col={'r.','rx','r*','rs','r^','rv'};

%plot numerical solution for different timestep sizes
for i = 1:numel(tau_range)
    
    %field_name for accessing the right values in struct
    field_name=sprintf('tau%i',i);
    
    %plot numerical solution and create legend entry
    h(i+1)=plot(numerical_solutions.(method_name).(field_name).t,...
        numerical_solutions.(method_name).(field_name).y);    
    leg_entry{i+1}=['numerical solution for tau=',...
        mat2str(numerical_solutions.(method_name).(field_name).tau)];
    
    %set color of plot, color is changing from [0 1 0] to [1 0 0] with
    %decreasing timestep size
    c=[(i-1)/(numel(tau_range)-1) 1-(i-1)/(numel(tau_range)-1) 0];
    set(h(i+1),'Color',c);
    
end

legend(h,leg_entry);
xlabel('t')
ylabel('p(t)')
axis([0 5 0 20]);
hold off


end

