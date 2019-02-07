function plotTotal(zone,tspan,x1,Uout,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD)
t1=tspan;
%zone=1;
lw=1.5;
lwt=12;
%scrsz = get(groot,'ScreenSize');
figure('Position',[50 50 600 700])
if zone==Leng+1
    subplot(4,3,1)
    plot(t1, x1(:,Leng+1),'LineWidth',lw);
    set(gca,'FontSize',12)
    %  title('Water temperature in Storage tank','Interpreter','latex')
    ylabel('$T_{st}$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
else
    subplot(4,3,1)
    plot(t1, x1(:,zone),'LineWidth',lw);
    set(gca,'FontSize',12)
    % title('Air temperature in Zones','Interpreter','latex')
    ylabel(['$T_{z_{ ' int2str(zone) ' }}$'],'Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
end
if zone==Leng+1
    subplot(4,3,4)
    plot(tout, Uout(:,Leng+1),'LineWidth',lw);
    %  title('Normalized energy in heat pump','Interpreter','latex')
    ylabel('$u_{st}$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    set(gca,'FontSize',12)
    grid on
else
    subplot(4,3,4)
    PosU=Leng+2:2*(Leng)+1;
    plot(tout, Uout(:,zone),'LineWidth',lw);
    set(gca,'FontSize',12)
    %  title('Flow rate in fan-coil units','Interpreter','latex')
    ylabel(['$u_{ ' int2str(zone) ' }$'],'Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
end
if zone==Leng+1
    subplot(4,3,7)
    plot(t1, x1(:,2*(Leng)+2),'LineWidth',lw);
    set(gca,'FontSize',12)
    %  title('Estimation of $T_{st}$','Interpreter','latex')
    ylabel('$\hat{T}_{st}$','Interpreter','latex','Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
else
    subplot(4,3,7)
    PosE=Leng+2:2*(Leng)+1;
    plot(t1, x1(:,PosE(zone)),'LineWidth',lw);
    set(gca,'FontSize',12)
    %  title(['Estimation of $T_{z_{' int2str(zone) '}}$'],'Interpreter','latex')
    ylabel(['$\hat{T}_{z_{' int2str(zone) ' }}$'],'Interpreter','latex', 'Fontsize',20)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    grid on
end 
if zone==Leng+1
    subplot(4,3,10)
    arr=[abs(Eout(:,zone)) Ebarout(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 55])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,'$\varepsilon_{y}^{s}$',...
        '$\overline{\varepsilon}_{y}^{s}$')
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='horizontal';
    h.Location='NorthEast';
else
    subplot(4,3,10)
    arr=[abs(Eout(:,zone)) Ebarout(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 22])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',15)
    ylabel(aa(2),['${D}^{(' int2str(zone) ')}$'],'interpreter',...
        'latex','Fontsize',20,'color',[0 0.6 0]);
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,['$| \varepsilon_{y}^{( ' int2str(zone) ' )}|$'],...
        ['$\overline{\varepsilon}_{y}^{( ' int2str(zone) ' )}$'])
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='horizontal';
    h.Location='NorthEast';
    %h.Position=[50 50 60 68];
end

if zone==Leng+1
    subplot(4,3,11)
    arr=[abs(Eouta(:,zone)) Ebarouta(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Douta(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 55])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
    for j=1:length(Eouta(:,zone))
         if tout(j) < TD(zone)
            bb(1).XData(j)=NaN;
            bb(1).YData(j)=NaN;
            bb(2).XData(j)=NaN;
            bb(2).YData(j)=NaN;
            cc.YData(j)=NaN;
         end
    end
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,'$\varepsilon_{y_a}^{s}$',...
        '$\overline{\varepsilon}_{y_a}^{s}$')
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='vertical';
    h.Location='NorthWest';
else
    subplot(4,3,11)
    arr=[abs(Eouta(:,zone)) Ebarouta(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Douta(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 30])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
    for j=1:length(Eouta(:,zone))
         if tout(j) < TD(zone)
            bb(1).XData(j)=NaN;
            bb(1).YData(j)=NaN;
            bb(2).XData(j)=NaN;
            bb(2).YData(j)=NaN;
            cc.YData(j)=NaN;
         end
    end
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',15)
    ylabel(aa(2),['${D}_{a}^{(' int2str(zone) ')}$'],'interpreter',...
        'latex','Fontsize',20,'color',[0 0.6 0]);
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,['$| \varepsilon_{y_a}^{( ' int2str(zone) ' )}|$'],...
        ['$\overline{\varepsilon}_{y_a}^{( ' int2str(zone) ' )}$'])
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='vertical';
    h.Location='NorthWest';
    %h.Position=[50 50 60 68];
end


if zone==Leng+1
    subplot(4,3,12)
    arr=[abs(Eouts(:,zone)) Ebarouts(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Douts(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 30])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
    for j=1:length(Eouts(:,zone))
         if tout(j) < TD(zone)
            bb(1).XData(j)=NaN;
            bb(1).YData(j)=NaN;
            bb(2).XData(j)=NaN;
            bb(2).YData(j)=NaN;
            cc.YData(j)=NaN;
         end
    end
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,'$\varepsilon_{y_s}^{s}$',...
        '$\overline{\varepsilon}_{y_s}^{s}$')
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')    
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='vertical';
    h.Location='NorthWest';
else
    subplot(4,3,12)
    arr=[abs(Eouts(:,zone)) Ebarouts(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Douts(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 30])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
%     for j=1:length(Eouts(:,zone))
%          if tout(j) < TD(zone)
%             bb(1).XData(j)=NaN;
%             bb(1).YData(j)=NaN;
%             bb(2).XData(j)=NaN;
%             bb(2).YData(j)=NaN;
%             cc.YData(j)=NaN;
%          end
%     end
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',15)
    ylabel(aa(2),['${D}_{s}^{(' int2str(zone) ')}$'],'interpreter',...
        'latex','Fontsize',20,'color',[0 0.6 0]);
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,['$| \varepsilon_{y_s}^{( ' int2str(zone) ' )}|$'],...
        ['$\overline{\varepsilon}_{y_s}^{( ' int2str(zone) ' )}$'])
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',20,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=20;
    h.Orientation='vertical';
    h.Location='NorthWest';
    %h.Position=[50 50 60 68];
end

subplot(4,3,2)
PosEa=2*Leng+3:3*Leng+2;
hh=plot(t1, x1(:,PosEa(zone)),'LineWidth',lw);
%     for j=1:length(x1(:,PosEa(zone)))
%          if t1(j) < TD(zone)
% 
%             hh.XData(j)=NaN;
%             hh.YData(j)=NaN;
%          end
%     end
set(gca,'FontSize',12)
ylabel(['$\hat{x}_{a}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
subplot(4,3,5)
PosEoa=3*Leng+3:4*Leng+2;
hh=plot(t1, x1(:,PosEoa(zone)),'LineWidth',lw);
    for j=1:length(x1(:,PosEoa(zone)))
         if t1(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
         end
    end
set(gca,'FontSize',12)
ylabel(['${\Omega}_{a}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
subplot(4,3,8)
PosEtha=4*Leng+3:5*(Leng)+2;
hh=plot(t1, x1(:,PosEtha(zone)),'LineWidth',lw);
    for j=1:length(x1(:,PosEtha(zone)))
         if t1(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
         end
    end
set(gca,'FontSize',12)
ylabel(['$\hat{f}_{a}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on

subplot(4,3,3)
PosEs=5*Leng+3:6*Leng+2;
hh=plot(t1, x1(:,PosEs(zone)),'LineWidth',lw);
    for j=1:length(x1(:,PosEs(zone)))
         if t1(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
         end
    end
set(gca,'FontSize',12)
ylabel(['$\hat{x}_{s}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
subplot(4,3,6)
PosEos=6*Leng+3:7*Leng+2;
hh=plot(t1, x1(:,PosEos(zone)),'LineWidth',lw);
    for j=1:length(x1(:,PosEos(zone)))
         if t1(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
         end
    end
set(gca,'FontSize',12)
ylabel(['${\Omega}_{s}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
subplot(4,3,9)
PosEths=7*Leng+3:8*(Leng)+2;
hh=plot(t1, x1(:,PosEths(zone)),'LineWidth',lw);
    for j=1:length(x1(:,PosEths(zone)))
         if t1(j) < TD(zone)
            hh.XData(j)=NaN;
            hh.YData(j)=NaN;
         end
    end
set(gca,'FontSize',12)
ylabel(['$\hat{f}_{s}^{(' int2str(zone) ' )}$'],'Interpreter','latex', 'Fontsize',20)
xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
grid on
end