function plotFault(zone,tspan,x1,tout,fout,fhatout,Df,Douts,Kout,Leng,TD)
t1=tspan;
%zone=1;
lw=2.5;
lwt=22;
mm=25;

figure('Position',[50 50 510 420])
%mm=plot(tout,fout(:,zone),'LineWidth',lw);
mm=plot(tout,Df(zone)*ones(1,length(tout)),'LineWidth',lw,'LineStyle','--');
hold on
PosEths=7*Leng+3:8*(Leng)+2;

a=abs( fout(:,zone)-fhatout(:,zone));

hh=plot(tout,a,'LineWidth',lw);
    for j=1:length(a)%x1(:,PosEths(zone))
         if tout(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
            mm.XData(j)=NaN;
            mm.YData(j)=NaN;
         end
    end
    for i=1:length(fout(:,zone))
         if tout(i) < TD(zone)-0.05
            mm.XData(i)=NaN;
            mm.YData(i)=NaN;
         end
    end
    for k=1:length(Douts(:,zone))
        if Douts(k,zone)==1
            Tisol=tout(k);
            break
        end
    end
    for g=1:length(Kout(:,zone))
        if Kout(g,zone)==1
            Tacco=tout(g);
            break
        end
    end
        
%set(gca,'FontSize',12)
ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15) 
h=legend(['$\Delta f^{(' int2str(zone) ' )}$'],['$|f^{(' int2str(zone) ' )}-\widehat{f}^{(' int2str(zone) ' )}|$']);%,'Interpreter','latex', 'Fontsize',20)
    h.Interpreter='latex';
    h.FontSize=30;
    h.Orientation='vertical';
    h.Location='NorthWest';
    h.Position=[0.240 0.487 0.310 0.365];%h.Position=[0.716 0.547 0.078 0.177];
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
xlim([0 t1(end)])
ylim([0 10])
set(gca,'XTick',[0 Tacco t1(end)],'FontSize',25,'LineWidth',1.5)%TD(zone)
title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
tightfig