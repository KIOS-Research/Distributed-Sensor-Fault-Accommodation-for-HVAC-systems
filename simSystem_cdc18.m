clc;
clear all
%close all
format long g
start_toolkit

% open this file to change the parameters of the simulations
collect_data_building_cdc18;

Leng=length(str_z);
a=1:(Leng);
xi(a)= str_z(a).Tz;


xi(Leng+1)=Twi;
xi(Leng+2:2*(Leng)+2)=0; % estimator length DETECTION
% ADAPTIVE ESTIMATION SCHEME FOR ACTUATOR FAULT IDENTIFICATION
for i=1:Leng
xi(2*Leng+2+i)=str_z(i).Tref; % adaptive estimation length
end
xi(3*(Leng)+3)=Trefw;
xi(3*Leng+4:4*(Leng)+4)=0; % omega length
xi(4*Leng+5:5*(Leng)+4)=0; % adaptation 
% ADAPTIVE ESTIMATION SCHEME FOR SENSOR FAULT IDENTIFICATION
xi(5*Leng+5:6*Leng+4)=0; % adaptive estimation length
xi(6*Leng+5:7*Leng+4)=0; % omega length
xi(7*Leng+5:8*(Leng)+4)=0; % adaptation 
dd(1:Leng+1)=0;



% initialization of the matrices
k=1:Leng;
for i=1:Leng
    for j=1:length(str_z(i).connected)
        con(i,j)=str_z(i).connected(j);
        Az(i,j)=str_z(i).az_ij(j);
    end
end
Umax(k)=str_z(k).Umax;
cz(k)=str_z(k).cz;
az(k)=str_z(k).az;
Ad(k)=str_z(k).Ad;
T1(k)=str_z(k).T1;

%initialization of the time series
Uout=zeros(0,Leng+1);
yout=zeros(0,Leng+1);
yVout=zeros(0,Leng+1);
fout=zeros(0,Leng+1);
fhatout=zeros(0,Leng+1);
faout=zeros(0,Leng+1);
Eout=zeros(0,Leng+1);
Ebarout=zeros(0,Leng+1);
Ebarouta=zeros(0,Leng);
Ebarouts=zeros(0,Leng);
Eouta=zeros(0,Leng);
Eouts=zeros(0,Leng+1);
Dout=zeros(0,Leng+1);
Kout=zeros(0,Leng);
Douta=zeros(0,Leng);
Douts=zeros(0,Leng);
Tout=zeros(0,0);
rout=zeros(0,Leng+1);
zout=zeros(0,Leng+1);
sout=zeros(0,Leng+1);
ebarout1=zeros(0,Leng);
ebarout2=zeros(0,Leng);
eFout=zeros(0,Leng);
guout=zeros(0,Leng+1);
gubarout=zeros(0,Leng);

Dgout=zeros(0,0);
gbarout=zeros(0,0);


ymax=40;
for i=1:(Leng)
    sum4=0;
    sum5=0;
    sum6=0;
    for c=1:(length(str_z(i).connected))     %sum for A
        sum4=sum4+str_z(i).az_ij(c)/str_z(i).cz;
        sum5=sum5+(str_z(i).az_ij(c)/str_z(i).cz)*...
            str_z(str_z(i).connected(c)).n_bar;  
        for d=1:str_z(i).paths.TotalPaths
            sum6=sum6+ ((p_air * Cp) / str_z(i).cz)*str_z(i).paths.Ad_ij(d)...
                *2*(ymax+max(str_z(i).n_bar,str_z(str_z(i).paths.ConnDoors(d)).n_bar))...
                *sqrt(str_z(i).n_bar+str_z(str_z(i).paths.ConnDoors(d)).n_bar);   
        end
    end
    A(i,i)=(h*str_z(i).Ad-(str_z(i).az)) /str_z(i).cz - sum4;
    
    va(i)=str_z(i).r_bar + abs(str_z(i).k)*str_z(i).n_bar...
        + (((str_z(i).Umax * awz)) /(str_z(i).cz))*(n_bar_s+str_z(i).n_bar)...
        + sum5 + sum6;
    
     Df(i)= abs( (va(i) + (A(i,i)-str_z(i).k)*str_z(i).n_bar) /...
         ( (((str_z(i).Umax * awz)) /(str_z(i).cz)) + A(i,i) ) );
     
     omega_ss(i)= abs( ( (((str_z(i).Umax * awz)) /(str_z(i).cz)) + abs(str_z(i).L) )/...
         ( A(i,i) - str_z(i).L) )
end
     
Df;
omega_ss;

% time parameters
tt=0;
step=0.01; % 0.02 %5e-3; %0.01; %2e-3; %8e-3;
time=1.5;
tspan=0:step:time; 

DET=1;
global TD TDa TDs TLs H PP KK W
W=10;
TD=(time+1)*ones(1,Leng+1);
TDa=(time+1)*ones(1,Leng);
TDs=(time+1)*ones(1,Leng);
TLs=(time+1)*ones(1,Leng);
H=zeros(W,Leng);
PP=zeros(1,Leng);
KK=zeros(1,Leng);
nn=8;

%% systemmv9 contains the differential equations
% need to be solve in order (1) to simulate the HVAC control system
% and also the diagnosis algorithm


tic
Tmin=ans

F=@(t,x) systemmv9(t, x, time, step, tspan, tt, nn, str_z, k, con,...
    Az, Umax, cz ,az, Ad, T1, h, awz, Cw, Ustmax, COPmax, To,...
    DTmax, Cp, Cv, p_air, Ta, Tpl, aw, maxst, minst, ...
    kst, Trefw, L_s, n_bar_s, p_s, lambda_s, x_bar_s,...
    r_bar_s, K, Leng, amp_out, amp_s, f_out, f_s, F_time_s, F_value_s,DET,Df);

x1= ode3(F,tspan, xi');
tout=Tout';

toc
Tmax=ans

