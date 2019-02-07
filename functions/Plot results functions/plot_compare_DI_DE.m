close all
clear all
clc
load WithoutAccommodation_z2_v1
zone=2;
plotTemp(zone,tspan,str_z,x1,Uout,tout,yout,Leng)
title(['$\Sigma^{(' int2str(zone) ')}$, without FA'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
%plotTotal(zone,tspan,x1,Uout,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD)    
%plotInput(zone,tspan,x1,Uout,tout,yout,Leng,str_z,Trefw)
plotARRs(zone,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD,str_z)
    
clear all
close all
clc
load Fault_Accommodation_z2_v2
zone=2;
plotTempV(zone,tspan,str_z,x1,Uout,tout,yout,yVout,Leng)
title(['$\Sigma^{(' int2str(zone) ')}$, FA'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
%plotTotal(zone,tspan,x1,Uout,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD)    
plotARRs(zone,tout,Eout,Ebarout,Dout,Kout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD,str_z)
plotFault(zone,tspan,x1,tout,fout,Dout,Douts,Kout,Leng,TD)

    
%%

load Distributed_control_8
%subplot(2,2,1)
plotARRs(zone,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD,str_z)
title(['$\Sigma^{(' int2str(zone) ')}$, DI'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')
load Decetralized_control_16
%subplot(2,2,2)
plotARRs(zone,tout,Eout,Ebarout,Dout,Eouta,Ebarouta,Douta,Eouts,Ebarouts,Douts,Leng,TD,str_z)
title(['$\Sigma^{(' int2str(zone) ')}$, DE'],'interpreter','latex',...
        'FontSize',27,'FontName','Times New Roman')