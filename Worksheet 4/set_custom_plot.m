function [] = set_custom_plot( vector_N_x, vector_tau, interesting_time, subplot_id )
%SET_CUSTOM_PLOT 

FONTSIZE=12;
for j = 1:numel(interesting_time)
    figure(j)
    set(gcf,'Position',[100 100 1120 520]);
    subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
    hold on
    set(gcf,'Color','white');
    text(0.1,0,['$t = ',mat2str(interesting_time(j)),'$'],'interpreter','latex','Fontsize',FONTSIZE);
    set(gca,'Color','white');
    set(gca,'XColor','white');
    set(gca,'YColor','white');
    hold off
end

subplot_id = 1;
for i = 1:numel(vector_tau)
    subplot_id = subplot_id +1;
    for j = 1:numel(interesting_time)
        figure(j)
        subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
        hold on
        set(gcf,'Color','white');       
        text(0.2,0,['$\tau = \frac{1}{',mat2str(1./vector_tau(i)),'}$'],'interpreter','latex','Fontsize',FONTSIZE);        
        set(gca,'Color','white');
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        hold off
    end
end

subplot_id=1;
for i = 1:numel(vector_N_x)
    subplot_id=subplot_id+1+numel(vector_tau);
    for j = 1:numel(interesting_time)
        figure(j)
        subplot(numel(vector_N_x)+1,numel(vector_tau)+1,subplot_id)
        hold on
        set(gcf,'Color','white');
        text(0.1,0.5,['$N_x = ',mat2str(vector_N_x(i)),'$'],'interpreter','latex','Fontsize',FONTSIZE);
        set(gca,'Color','white');
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        hold off
    end
end


end

