%load cdc18_v11
for zone=1:3
    
    
x_size=32;
tick_size=25;
y_size=35;
t_size=40;

if zone==84
figure
hold on
plot(tout, yout(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tspan, x1(:,zone),'LineWidth',2.5,'color','[0.7 0 0]');
plot(tspan, x1(:,2*Leng+2+zone),'LineWidth',2.5,'color','[0 0.6 0.86]'); 
plot(tspan, Trefw*ones(1,length(tspan)),'LineWidth',2.5,'color','[0.3 0.3 0.3]','LineStyle','--');
set(gca,'FontSize',t_size)
h=legend('$y^{s}$','$T_{st}$','$\widehat{T}_{st}$',...
     '$y^{s}_{r}$');
    h.Interpreter='latex';
    h.FontSize=t_size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',t_size)
    grid on
    fname=['Temp_' int2str(zone)];
    saveas(gcf,fname,'pdf');

 figure
hold on
plot(tout, Eouts(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tout, fout(:,zone),'LineWidth',2.5,'color','[0.7 0 0]')
%plot(tout, fhatout(:,zone),'LineWidth',2.5,'LineStyle','--')
plot(tspan, x1(:,3*Leng+3+zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',t_size)
h=legend('$\varepsilon^{s}$',...
    '$f^{s}$','$\widehat{f}^{s}$');
    h.Interpreter='latex';
    h.FontSize=t_size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',t_size)   
    grid on
    fname=['fault_' int2str(zone)];
    saveas(gcf,fname,'pdf');  
    
figure
plot(tout, Uout(:,zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',t_size)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$u_{st}$','Interpreter','latex','Fontsize',y_size)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',t_size)  
    grid on
    fname=['u_' int2str(zone)];
    saveas(gcf,fname,'pdf'); 
    
figure
plot(tout, yout(:,zone)-str_z(zone).Tref*ones(1,length(tout)),...
    'LineWidth',2.5,'color','[0.5 0.9 0.86]')%,'LineStyle',':');
set(gca,'FontSize',t_size)
hold on
t=x1(:,zone)-str_z(zone).Tref*ones(1,length(tspan));
plot(tspan, t,'LineWidth',2.5,'color','[0.7 0 0]');
h=legend(['$y^{(' int2str(zone) ')}-y^{(' int2str(zone) ')}_{r}$'],...
    ['$\epsilon^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=t_size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',t_size) 
    grid on    
    fname=['e_' int2str(zone)];
    saveas(gcf,fname,'pdf'); 
    
else

figure
hold on
plot(tout, zout(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tspan, x1(:,zone),'LineWidth',2.5,'color','[0.7 0 0]');
%plot(tspan, x1(:,2*Leng+2+zone),'LineWidth',2.5,'color','[0 0.6 0.86]'); 
plot(tspan, str_z(zone).Tref*ones(1,length(tspan)),'LineWidth',2.5,'color','[0.3 0.3 0.3]','LineStyle','--');
set(gca,'XTick',0:0.2:1.5,'FontSize',tick_size)
h=legend(['$z^{(' int2str(zone) ')}$'],['$T_{z_' int2str(zone) '}$'],...
     ['$y^{(' int2str(zone) ')}_{r}$']);
    h.Interpreter='latex';
    h.FontSize=t_size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',t_size)
    grid on
    fname=['Temp_' int2str(zone)];
    saveas(gcf,fname,'pdf');
    
    
 figure
hold on
%plot(tout, Eouts(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tout, fout(:,zone),'LineWidth',2.5,'color','[0.7 0 0]')
%plot(tout, fhatout(:,zone),'LineWidth',2.5,'LineStyle','--')
plot(tspan, x1(:,3*Leng+3+zone),'LineWidth',2.5,'color','[0 0.6 0.86]','LineStyle','--');
set(gca,'XTick',0:0.2:1.5,'FontSize',tick_size)
h=legend(['$f^{(' int2str(zone) ')}$'],['$\widehat{f}^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=t_size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',t_size)   
    grid on
    fname=['fault_' int2str(zone)];
    saveas(gcf,fname,'pdf');
   
figure
plot(tout, Uout(:,zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'XTick',0:0.2:1.5,'FontSize',tick_size)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    ylabel(['$u_{' int2str(zone) '}$'],'Interpreter','latex','Fontsize',y_size)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',t_size)  
    grid on
    fname=['u_' int2str(zone)];
    saveas(gcf,fname,'pdf');
    
    
figure
p=yout(:,zone)-str_z(zone).Tref*ones(1,length(tout(1)));
% plot(tout,p ,...
%     'LineWidth',2.5,'color','[0.5 0.9 0.86]')%,'LineStyle',':');
plot(tout,sout(:,zone) ,...
    'LineWidth',2.5,'color','[0.5 0.9 0.86]')
hold on
k=x1(:,zone)-str_z(zone).Tref*ones(1,length(tspan(1)));
plot(tspan, k,'LineWidth',2.5,'color','[0.7 0 0]');
%plot(tout,sout(:,zone),'LineWidth',2.5,'color','[0 0.6 0.86]')

h=legend(['$\epsilon_{y}^{(' int2str(zone) ')}-\widehat{f}^{(' int2str(zone) ')}$'],...
    ['$\epsilon^{(' int2str(zone) ')}$'])
%,...
 %   ['$\epsilon_{y}^{(' int2str(zone) ')}-\widehat{f}^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',x_size)
    set(gca,'XTick',0:0.2:1.5,'FontSize',tick_size)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',y_size)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',t_size) 
    grid on
    fname=['e_' int2str(zone)];% 
    saveas(gcf,fname,'pdf'); 
end

end