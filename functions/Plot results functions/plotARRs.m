function plotARRs(zone,tout,Eout,Ebarout,Dout,Kout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD,str_z)
lwt=12;
lw=2;
lwD=2.2;
mm=25;

if zone==Leng+1
    figure('Position',[50 50 510 420]) %[50 50 520 480]
%     daspect([5.6 5.5 1])
    arr=[abs(Eout(:,zone)) Ebarout(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 30])
    ylim(aa(2),[0 1.02])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lwD)
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(aa(1),'$\varepsilon_{y}^{s}$',...
        '$\overline{\varepsilon}_{y}^{s}$','$D^{s}$')
    %title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
    %    'FontSize',mm,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='NorthWest';
    h.Position=[0.240 0.487 0.310 0.365];%h.Position=[0.716 0.547 0.078 0.177];
            set(aa(1),'YTick',0:10:80,'FontSize',mm,'LineWidth',1.2)
    set(aa(2),'YTick',0:1,'FontSize',mm)
        title(['$\mathcal{M}^{s}$'],'interpreter','latex',...
        'FontSize',25,'FontName','Times New Roman')

    tightfig
    fname='arr_s';
    saveas(gcf,fname,'eps2c');
else
    figure('Position',[50 50 510 420])
    arr=[abs(Eout(:,zone)) Ebarout(:,zone)];
    %Decision=[Dout(:,zone) Kout(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
    %[aa,bb,cc]=plotyy(tout,arr,tout,Decision);
    %hobj=plotyy(tout,arr,tout,Dout(:,zone));
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 8])
    ylim(aa(2),[0 1.02])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lwD)
    set(aa(1),'YTick',0:2:22,'FontSize',35,'LineWidth',1.2)
    set(aa(2),'YTick',0:1,'FontSize',35)
    doit=0;
    for i=1:length(cc.YData)
        if doit==0 && cc.YData(i)==1
            yclick=cc.YData(i)
            xclick=cc.XData(i)
            doit=1;
            break
        end
    end
    %datacursormode on
 %   makedatatip(cc,[i i])
    if doit==1
        hDatatip = zeros(1,1);
        hDataCursorMgr =datacursormode(ancestor(cc,'figure'))
        hDatatip = createDatatip(hDataCursorMgr, cc);
        pos = [xclick yclick];
%         set(get(hDatatip),'DataIndex',i,...
%               'TargetPoint',pos)  
        set(hDatatip,'Position',pos)
        info_struct = getCursorInfo(hDataCursorMgr)
        set(info_struct.Target,'LineWidth',2.5) 
        updateDataCursors(hDataCursorMgr)
        alldatacursors = findall(gcf,'type','hggroup')
        set(alldatacursors,'FontSize',25)
        set(alldatacursors,'FontName','Times')
    end

    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
%     set(aa(1),'dataaspectratio',[2 4 1])
%     set(aa(2),'dataaspectratio',[1 1 1])
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
    ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',mm)
    %ylabel(aa(2),['${D}^{(' int2str(zone) ')}$'],'interpreter',...
    %    'latex','Fontsize',20,'color',[0 0.6 0]);
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(aa(1),['$| \varepsilon_{y}^{( ' int2str(zone) ' )}|$'],...
        ['$\overline{\varepsilon}_{y}^{( ' int2str(zone) ' )}$'],...
        ['${D}^{(' int2str(zone) ')}$'])
    h.Interpreter='latex';
    h.FontSize=30;
    h.Orientation='vertical';
    h.Location='NorthWest';
    h.Position=[0.233 0.474 0.318 0.365];
            set(aa(1),'YTick',0:2:22,'FontSize',mm,'LineWidth',1.2)
    set(aa(2),'YTick',0:1,'FontSize',mm)
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
    tightfig
%     fname=['arr_' int2str(zone)];
%     saveas(gcf,fname,'eps2c');
end


% if zone==Leng+1
%     figure('Position',[50 50 520 500])
%     arr=[abs(Eouta(:,zone)) Ebarouta(:,zone)];
%     [aa,bb,cc]=plotyy(tout,arr,tout,Douta(:,zone));
%     set(aa,{'ycolor'},{'k';[0 0.6 0]})
%     ylim(aa(1),[-0.05 55])
%     ylim(aa(2),[0 1.05])
%     grid on
%     set(bb(1),'LineWidth',lw)
%     set(bb(2),'LineWidth',lw,'LineStyle','--')
%     set(cc,'LineWidth',lw,'color',[0 0.6 0])
%     set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lwD)
%     for j=1:length(Eouta(:,zone))
%          if tout(j) < TD(zone)
%             bb(1).XData(j)=NaN;
%             bb(1).YData(j)=NaN;
%             bb(2).XData(j)=NaN;
%             bb(2).YData(j)=NaN;
%             cc.YData(j)=NaN;
%          end
%     end
%     %     alldatacursors = findall(gcf,'type','hggroup')
%     %     set(alldatacursors,'FontSize',20)
%     % axis font size
%     set(aa(1),'FontSize',12)
%     set(aa(2),'FontSize',12)
%     set(aa(1),'YTick',0:2:22,'FontSize',35,'LineWidth',1.2)
%     set(aa(2),'YTick',0:1,'FontSize',35)
%     ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
%     xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
%     h=legend(aa(1),'$\varepsilon_{y_a}^{s}$',...
%         '$\overline{\varepsilon}_{y_a}^{s}$')
%     title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
%         'FontSize',mm,'FontName','Times New Roman')
%     h.Interpreter='latex';
%     h.FontSize=21;
%     h.Orientation='vertical';
%     h.Location='NorthWest';
%     set(gca,'fontsize',mm)
% else
%     figure('Position',[50 50 520 500])
%     arr=[abs(Eouta(:,zone)) Ebarouta(:,zone)];
%     [aa,bb,cc]=plotyy(tout,arr,tout,Douta(:,zone));
%     set(aa,{'ycolor'},{'k';[0 0.6 0]})
%     ylim(aa(1),[-0.05 22])
%     ylim(aa(2),[0 1.05])
%     grid on
%     set(bb(1),'LineWidth',lw)
%     set(bb(2),'LineWidth',lw,'LineStyle','--')
%     set(cc,'LineWidth',lw,'color',[0 0.6 0])
%     set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lwD)
%     for j=1:length(Eouta(:,zone))
%          if tout(j) < TD(zone)
%             bb(1).XData(j)=NaN;
%             bb(1).YData(j)=NaN;
%             bb(2).XData(j)=NaN;
%             bb(2).YData(j)=NaN;
%             cc.YData(j)=NaN;
%          end
%     end
%     %     alldatacursors = findall(gcf,'type','hggroup')
%     %     set(alldatacursors,'FontSize',20)
%     % axis font size
%         set(aa(1),'YTick',0:2:22,'FontSize',35,'LineWidth',1.2)
%     set(aa(2),'YTick',0:1,'FontSize',35)
%     set(aa(1),'FontSize',12)
%     set(aa(2),'FontSize',12)
%     ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',15)
%     %ylabel(aa(2),['${D}_{a}^{(' int2str(zone) ')}$'],'interpreter',...
%     %    'latex','Fontsize',20,'color',[0 0.6 0]);
%     xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
%     h=legend(aa(1),['$| \varepsilon_{y_a}^{( ' int2str(zone) ' )}|$'],...
%         ['$\overline{\varepsilon}_{y_a}^{( ' int2str(zone) ' )}$'],...
%         ['${I}_{a}^{(' int2str(zone) ')}$'])
%     title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
%         'FontSize',22,'FontName','Times New Roman')
%         set(aa(1),'YTick',0:2:22,'FontSize',mm,'LineWidth',1.2)
%     set(aa(2),'YTick',0:1,'FontSize',mm)
%     h.Interpreter='latex';
%     h.FontSize=21;
%     h.Orientation='vertical';
%     h.Location='NorthWest';
%     %h.Position=[50 50 60 68];
%     set(gca,'fontsize',mm)
% end
% 
% 
if zone==Leng+1
    figure('Position',[50 50 510 420])
    arr=[abs(Eouts(:,zone)) Ebarouts(:,zone)];
    %[aa,bb,cc]=plotyy(tout,arr,tout,Douts(:,zone));
    Decision=[Douts(:,zone) Kout(:,zone)];
    %[aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
    [aa,bb,cc]=plotyy(tout,arr,tout,Decision);
    set(aa,{'ycolor'},{'k';[0 0.6 0]})
    ylim(aa(1),[-0.05 22])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    set(cc,'LineWidth',lw,'color',[0 0.6 0])
    set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lwD)
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
        set(aa(1),'YTick',0:2:22,'FontSize',35,'LineWidth',1.2)
    set(aa(2),'YTick',0:1,'FontSize',35)
    ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(bb,'$\varepsilon_{y_s}^{s}$',...
        '$\overline{\varepsilon}_{y_s}^{s}$')
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',mm,'FontName','Times New Roman')    
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='NorthWest';
    set(gca,'fontsize',mm)
else
    figure('Position',[50 50 510 420])
    arr=[abs(Eouts(:,zone)) Ebarouts(:,zone)];
    Decision=[Douts(:,zone) Kout(:,zone)];
    [aa,bb,cc]=plotyy(tout,arr,tout,Decision);
    set(aa,{'ycolor'},{'k';[0 0 0]})
    ylim(aa(1),[-0.05 8])
    ylim(aa(2),[0 1.05])
    grid on
    set(bb(1),'LineWidth',lw)
    set(bb(2),'LineWidth',lw,'LineStyle','--')
    %set(cc(1),'LineWidth',lw,'color',[0 0.6 0])
    set(cc(1),'color',[0 0.6 0],'LineStyle','-.','LineWidth',lwD)
    %set(cc(2),'LineWidth',lw,'color',[0 1 0])
    set(cc(2),'color',[0.5 0 0.5],'LineStyle',':','LineWidth',lwD)
    for j=1:length(Eouts(:,zone))
         if tout(j) < TD(zone)
            bb(1).XData(j)=NaN;
            bb(1).YData(j)=NaN;
            bb(2).XData(j)=NaN;
            bb(2).YData(j)=NaN;
            cc(1).YData(j)=NaN;
            cc(2).YData(j)=NaN;
         end
    end
    doit=0;
    for i=1:length(cc(2).YData)
        if doit==0 && cc(2).YData(i)==1
            yclick=cc(2).YData(i)
            xclick=cc(2).XData(i)
            doit=1;
            break
        end
    end
    if doit==1
        hDatatip = zeros(1,1);
        hDataCursorMgr =datacursormode(ancestor(cc(2),'figure'))
        hDatatip = createDatatip(hDataCursorMgr, cc(2));
        pos = [xclick yclick];
%         set(get(hDatatip),'DataIndex',i,...
%               'TargetPoint',pos)  
        set(hDatatip,'Position',pos)
        info_struct = getCursorInfo(hDataCursorMgr)
        set(info_struct.Target,'LineWidth',2.5) 
        updateDataCursors(hDataCursorMgr)
        alldatacursors = findall(gcf,'type','hggroup')
        set(alldatacursors,'FontSize',25)
        set(alldatacursors,'FontName','Times')
    end
    
    %     alldatacursors = findall(gcf,'type','hggroup')
    %     set(alldatacursors,'FontSize',20)
    % axis font size
    set(aa(1),'FontSize',12)
    set(aa(2),'FontSize',12)
        set(aa(1),'YTick',0:2:22,'FontSize',mm,'LineWidth',1.2)
    set(aa(2),'YTick',0:1,'FontSize',mm)
    ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',15)
    %ylabel(aa(2),['${D}_{s}^{(' int2str(zone) ')}$'],'interpreter',...
    %    'latex','Fontsize',20,'color',[0 0.6 0]);
    xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
    h=legend(aa(1),['$| \varepsilon_{y_s}^{( ' int2str(zone) ' )}|$'],...
        ['$\overline{\varepsilon}_{y_s}^{( ' int2str(zone) ' )}$'],...
        ['${I}_{s}^{(' int2str(zone) ')}$'],['${A}_{s}^{(' int2str(zone) ')}$'])
    title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
        'FontSize',mm,'FontName','Times New Roman')
    h.Interpreter='latex';
    h.FontSize=30;
    h.Orientation='vertical';
    h.Location='NorthWest';
    %h.Position=[50 50 60 68];
    set(gca,'fontsize',mm)
end
tightfig

%% plot agents in one figure

% mm=18;
% marker=['o','d','*','^'];
% figure('Position',[50 50 520 500])
% hold all
% for zone=2:5
% if zone==Leng+1
%     arr=[abs(Eout(:,zone)) Ebarout(:,zone)];
%     [aa,bb,cc]=plotyy(tout,arr,tout,Dout(:,zone));
%     set(aa,{'ycolor'},{'k';[0 0.6 0]})
%     ylim(aa(1),[-0.05 55])
%     ylim(aa(2),[0 1.05])
%     grid on
%     set(bb(1),'LineWidth',lw)
%     set(bb(2),'LineWidth',lw,'LineStyle','--')
%     set(cc,'LineWidth',lw,'color',[0 0.6 0])
%     set(cc,'color',[0 0.6 0],'LineStyle',':','LineWidth',lw)
%     %     alldatacursors = findall(gcf,'type','hggroup')
%     %     set(alldatacursors,'FontSize',20)
%     % axis font size
%     set(aa(1),'FontSize',12)
%     set(aa(2),'FontSize',12)
%     ylabel('$(^\circ C)$','Interpreter','latex','Fontsize',15)
%     xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
%     h=legend(bb,'$\varepsilon_{y}^{s}$',...
%         '$\overline{\varepsilon}_{y}^{s}$')
%     title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
%         'FontSize',mm,'FontName','Times New Roman')
%     h.Interpreter='latex';
%     h.FontSize=20;
%     h.Orientation='vertical';
%     h.Location='NorthWest';
%     set(gca,'fontsize',mm)
% else
%     arr=[abs(Eout(1:40+zone:end,zone)) Ebarout(1:40+zone:end,zone)];
%     [aa,bb,cc]=plotyy(tout(1:40+zone:end),arr,tout(1:40:end),Dout(1:40:end,zone));
%     set(aa,{'ycolor'},{'k';[0 0.6 0]})
%     ylim(aa(1),[-0.05 22])
%     ylim(aa(2),[0 1.05])
%     grid on
%     
%     set(bb(1),'LineWidth',lw,'Marker',marker(zone-1),'color',[0 0.447 0.741])
%     set(bb(2),'LineWidth',lw,'Marker',marker(zone-1),'color',[0.85 0.325 0.098])
%     %set(cc(1:30:end),'LineWidth',lw,'color',[0 0.6 0])
%     set(cc(1:30:end),'color',[0 0.6 0],'LineWidth',lw,'Marker',marker(zone-1))
%     set(aa(1),'YTick',0:2:22,'FontSize',35,'LineWidth',1.2)
%     set(aa(2),'YTick',0:1,'FontSize',35)
%     %     alldatacursors = findall(gcf,'type','hggroup')
%     %     set(alldatacursors,'FontSize',20)
%     % axis font size
%     set(aa(1),'FontSize',12)
%     set(aa(2),'FontSize',12)
%     ylabel(aa(1),'$(^\circ C)$','Interpreter','latex','Fontsize',mm)
%     %ylabel(aa(2),['${D}^{(' int2str(zone) ')}$'],'interpreter',...
%     %    'latex','Fontsize',20,'color',[0 0.6 0]);
%     xlabel('Time (hours)','Interpreter','latex','Fontsize',lwt)
% %     h=legend(aa(1),['$| \varepsilon_{y}^{( ' int2str(zone) ' )}|$'],...
% %         ['$\overline{\varepsilon}_{y}^{( ' int2str(zone) ' )}$'],...
% %         ['${D}^{(' int2str(zone) ')}$'])
% %     h.Interpreter='latex';
% %     h.FontSize=21;
% %     h.Orientation='vertical';
% %     h.Location='NorthWest';
%     %h.Position=[50 50 60 68];
%     title(['$\mathcal{M}^{(' int2str(zone) ')}$'],'interpreter','latex',...
%         'FontSize',22,'FontName','Times New Roman')
%         set(aa(1),'YTick',0:2:22,'FontSize',mm,'LineWidth',1.2)
%     set(aa(2),'YTick',0:1,'FontSize',mm)
%     set(gca,'fontsize',mm)
% end
% end
% 
%     h=legend('$| \varepsilon_{y}^{(2)}|$','$\overline{\varepsilon}_{y}^{(2)}$','${D}^{(2)}$',...
%         '$| \varepsilon_{y}^{(3)}|$','$\overline{\varepsilon}_{y}^{(3)}$','${D}^{(3)}$',...
%         '$| \varepsilon_{y}^{(4)}|$','$\overline{\varepsilon}_{y}^{(4)}$','${D}^{(4)}$',...
%         '$| \varepsilon_{y}^{(5)}|$','$\overline{\varepsilon}_{y}^{(5)}$','${D}^{(5)}$')
%     h.Interpreter='latex';
%     h.FontSize=20;
%     h.Orientation='vertical';
%     h.Location='northeast';
end