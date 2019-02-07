function plotTempV(zone,tspan,str_z,x1,Uout,tout,yout,yVout,Leng)
t1=tspan;
%zone=1;
lw=2.2;
lwt=12;
mm=25;
%scrsz = get(groot,'ScreenSize');
figure('Position',[50 50 510 420])

if zone==Leng+1
    hold on
    plot(t1(1:2:length(t1)),55*ones(length(t1(1:2:length(t1))),1),'LineStyle','--','LineWidth',lw,'color','k')
    plot(tout,yout(:,84),'LineWidth',lw,'color','[0 0.7 0.5]')
    plot(t1, x1(:,Leng+1),'LineWidth',lw,'LineStyle','-.','color','[0.7 0 0]');
    plot(t1, x1(:,2*(Leng)+2),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]');
    plot(tout,yVout(:,84),'LineWidth',lw,'color','[0 0.7 0]')
    ylabel('$(^oC)$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
    
    h=legend('$y^{s}_{ref}$','$y^{s}$','$T_{st}$','$\hat{T}_{st}$')
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='SouthWest';
    h.Position=[0.373 0.264 0.071 0.246];
    set(gca,'YTick',0:5:80,'FontSize',mm,'LineWidth',1.2)
    title('$\Sigma^{s}$','interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
    
    ylim(gca,[30 60])
    tightfig
    fname='T_s';
    saveas(gcf,fname,'eps2c')
else
    hold on
    PosE=Leng+2:2*(Leng)+1;
    plot(t1(1:2:length(t1)),str_z(zone).Tref*ones(length(t1(1:2:length(t1))),1),'LineStyle','--','LineWidth',lw,'color','k')
    %plot(tout,yout(:,zone),'LineWidth',lw,'color','[0 0.7 0.7]')
    plot(tout,yVout(:,zone),'LineWidth',lw,'color','[0.7 0.9 0.9]')
    plot(t1, x1(:,zone),'LineWidth',lw,'color','[0.7 0 0]');
    plot(t1, x1(:,PosE(zone)),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]'); 
    PosEs=5*Leng+3:6*Leng+2;
    plot(t1, x1(:,PosEs(zone)),'LineWidth',lw,'color','[0 0.9 0]'); 
    
    
    ylabel('$(^oC)$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
    
    h=legend(['$y^{(' int2str(zone) ')}_{ref}$'],...
        ['$y^{(' int2str(zone) ')}_{v}$'],['$T_{z_{ ' int2str(zone) ' }}$'],['$\hat{T}_{z_{' int2str(zone) ' }}$'],['$\hat{T}^{I}_{z_{' int2str(zone) ' }}$'])
    h.Interpreter='latex';
    h.FontSize=27;
    h.Orientation='vertical';
    h.Location='SouthWest';
    h.Position=[0.373 0.264 0.071 0.246];
    set(gca,'YTick',0:1:80,'FontSize',mm,'LineWidth',1.2)
    title(['$\Sigma^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
    
    ylim(gca,[17 27])
    tightfig
%     fname=['T_' int2str(zone)];
%     saveas(gcf,fname,'eps2c')
end

end