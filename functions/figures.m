%load cdc18_v11
for zone=1:3
    
size=30;
if zone==84
figure
hold on
plot(tout, yout(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tspan, x1(:,zone),'LineWidth',2.5,'color','[0.7 0 0]');
plot(tspan, x1(:,2*Leng+2+zone),'LineWidth',2.5,'color','[0 0.6 0.86]'); 
plot(tspan, Trefw*ones(1,length(tspan)),'LineWidth',2.5,'color','[0.3 0.3 0.3]','LineStyle','--');
set(gca,'FontSize',18)
h=legend('$y^{s}$','$T_{st}$','$\widehat{T}_{st}$',...
     '$y^{s}_{r}$');
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',size)
    grid on
    fname=['Temp_' int2str(zone)];
    saveas(gcf,fname,'pdf');

 figure
hold on
plot(tout, Eouts(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tout, fout(:,zone),'LineWidth',2.5,'color','[0.7 0 0]')
%plot(tout, fhatout(:,zone),'LineWidth',2.5,'LineStyle','--')
plot(tspan, x1(:,3*Leng+3+zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',18)
h=legend('$\varepsilon^{s}$',...
    '$f^{s}$','$\widehat{f}^{s}$');
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',size)   
    grid on
    fname=['fault_' int2str(zone)];
    saveas(gcf,fname,'pdf');  
    
figure
plot(tout, Uout(:,zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',18)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$u_{st}$','Interpreter','latex','Fontsize',30)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',size)  
    grid on
    fname=['u_' int2str(zone)];
    saveas(gcf,fname,'pdf'); 
    
figure
plot(tout, yout(:,zone)-str_z(zone).Tref*ones(1,length(tout)),...
    'LineWidth',2.5,'color','[0.5 0.9 0.86]')%,'LineStyle',':');
set(gca,'FontSize',18)
hold on
t=x1(:,zone)-str_z(zone).Tref*ones(1,length(tspan));
plot(tspan, t,'LineWidth',2.5,'color','[0.7 0 0]');
h=legend(['$y^{(' int2str(zone) ')}-y^{(' int2str(zone) ')}_{r}$'],...
    ['$\epsilon^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title('$\Sigma^{s}$','Interpreter','latex','Fontsize',size) 
    grid on    
    fname=['e_' int2str(zone)];
    saveas(gcf,fname,'pdf'); 
    
else

figure
hold on
plot(tout, yout(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tspan, x1(:,zone),'LineWidth',2.5,'color','[0.7 0 0]');
plot(tspan, x1(:,2*Leng+2+zone),'LineWidth',2.5,'color','[0 0.6 0.86]'); 
plot(tspan, str_z(zone).Tref*ones(1,length(tspan)),'LineWidth',2.5,'color','[0.3 0.3 0.3]','LineStyle','--');
set(gca,'FontSize',18)
h=legend(['$y^{(' int2str(zone) ')}$'],['$T_{z_' int2str(zone) '}$'],['$\widehat{T}_{z_' int2str(zone) '}$'],...
     ['$y^{(' int2str(zone) ')}_{r}$']);
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',size)
    grid on
    fname=['Temp_' int2str(zone)];
    saveas(gcf,fname,'pdf');
    
    
 figure
hold on
plot(tout, Eouts(:,zone),'LineWidth',2.5,'color','[0.5 0.9 0.86]');
plot(tout, fout(:,zone),'LineWidth',2.5,'color','[0.7 0 0]')
%plot(tout, fhatout(:,zone),'LineWidth',2.5,'LineStyle','--')
plot(tspan, x1(:,3*Leng+3+zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',18)
h=legend(['$\varepsilon^{(' int2str(zone) ')}$'],...
    ['$f^{(' int2str(zone) ')}$'],['$\widehat{f}^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',size)   
    grid on
    fname=['fault_' int2str(zone)];
    saveas(gcf,fname,'pdf');
   
figure
plot(tout, Uout(:,zone),'LineWidth',2.5,'color','[0 0.6 0.86]');
set(gca,'FontSize',18)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    ylabel(['$u_{' int2str(zone) '}$'],'Interpreter','latex','Fontsize',30)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',size)  
    grid on
    fname=['u_' int2str(zone)];
    saveas(gcf,fname,'pdf');
    
    
figure
p=yout(:,zone)-str_z(zone).Tref*ones(1,length(tout(1)));
plot(tout,p ,...
    'LineWidth',2.5,'color','[0.5 0.9 0.86]')%,'LineStyle',':');
hold on
k=x1(:,zone)-str_z(zone).Tref*ones(1,length(tspan(1)));
plot(tspan, k,'LineWidth',2.5,'color','[0.7 0 0]');
h=legend(['$\epsilon_{y}^{(' int2str(zone) ')}$'],...
    ['$\epsilon^{(' int2str(zone) ')}$']);
    h.Interpreter='latex';
    h.FontSize=size;
    h.Orientation='vertical';
    h.Location='NorthWest';
    xlabel('Time (hours)','Interpreter','latex','Fontsize',20)
    set(gca,'FontSize',18)
    ylabel('$(^o$C)','Interpreter','latex','Fontsize',22)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'Interpreter','latex','Fontsize',size) 
    grid on
    fname=['e_' int2str(zone)];
    saveas(gcf,fname,'pdf'); 
end

end