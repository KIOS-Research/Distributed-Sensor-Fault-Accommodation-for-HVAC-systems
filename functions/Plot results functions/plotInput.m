function plotInput(zone,tspan,x1,Uout,tout,yout,Leng,str_z,Trefw)
t1=tspan;
%zone=1;
lw=2.2;
lwt=12;
mm=18;
%scrsz = get(groot,'ScreenSize');
figure('Position',[50 50 400 420])

if zone==Leng+1
    hold on
%     plot(t1(1:2:length(t1)),Trefw*ones(length(t1(1:2:length(t1))),1),'LineStyle','--','LineWidth',lw,'color','k')
%     plot(tout,yout(:,84),'LineWidth',lw,'color','[0 0.7 0.5]')
%     plot(t1, x1(:,Leng+1),'LineWidth',lw,'LineStyle','-.','color','[0.7 0 0]');
%     plot(t1, x1(:,2*(Leng)+2),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]');
  %  tracking_error=Trefw*ones(length(t1),1)-x1(:,Leng+1);
    plot(tout, Uout(:,Leng),'LineWidth',lw,'color','b')
    
    %ylabel('$(^oC)$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
    
    h=legend(['$u_{st}$'])
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='SouthWest';
    h.Position=[0.373 0.264 0.071 0.246];
    set(gca,'YTick',0:1:10,'FontSize',mm,'LineWidth',1.2)
    title('$\Sigma^{s}$','interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
    
    ylim(gca,[0 6])
    tightfig
    fname='T_s';
    saveas(gcf,fname,'eps2c')
else
    hold on
    PosE=Leng+2:2*(Leng)+1;
%     plot(t1(1:2:length(t1)),str_z(zone).Tref*ones(length(t1(1:2:length(t1))),1),'LineStyle','--','LineWidth',lw,'color','k')
%     plot(tout,yout(:,zone),'LineWidth',lw,'color','[0 0.7 0.5]')
%     plot(t1, x1(:,zone),'LineWidth',lw,'LineStyle','-.','color','[0.7 0 0]');
%     %plot(t1, x1(:,PosE(zone)),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]');
%    plot(t1(1:2:length(t1)),22.5*ones(length(t1(1:2:length(t1))),1),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]')
     
    %tracking_error= x1(:,zone)-str_z(zone).Tref*ones(length(t1),1);
    %output_tracking_error= yout(:,zone)-str_z(zone).Tref*ones(length(tout),1); 
    
    %plot(tout, output_tracking_error,'LineWidth',lw,'color','[0.5 0.5 0.5]')
    %plot(t1, tracking_error,'LineWidth',lw,'color','b')
    plot(tout, Uout(:,zone),'LineWidth',lw,'color','b')
%     plot(t1(1:2:length(t1)),(-str_z(zone).Tref+25)*ones(length(t1(1:2:length(t1))),1),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]')
%     plot(t1(1:2:length(t1)),(-str_z(zone).Tref+22.5)*ones(length(t1(1:2:length(t1))),1),'LineWidth',lw,'LineStyle',':','color','[0 0 0.7]')
    %ylabel('$(^oC)$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
    
%     h=legend(['$y^{(' int2str(zone) ')}_{ref}$'],['$y^{(' int2str(zone) ')}$'],...
%         ['$T_{z_{ ' int2str(zone) ' }}$'],['$\hat{T}_{z_{' int2str(zone) ' }}$'],['$[{T}^{\min}_{z_{' int2str(zone) ' }},{T}^{\max}_{z_{' int2str(zone) ' }}]$'])
    h=legend(['$u_{' int2str(zone) '}$'])
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='SouthWest';
    h.Position=[0.373 0.264 0.071 0.246];
    set(gca,'YTick',-1:0.1:2,'FontSize',mm,'LineWidth',1.2)
    title(['$\Sigma^{(' int2str(zone) ')}, K^{(' int2str(zone) ')}_{H}=' int2str(str_z(zone).k) '$'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
%     title(['$\Sigma^{(' int2str(zone) ')},$ FTC'],'interpreter','latex',...
%         'FontSize',27,'FontName','Times New Roman')
    ylim(gca,[-0.01 1.01])
    tightfig
    fname=['tildet_x_' int2str(zone)];
    saveas(gcf,fname,'eps2c')
end

end