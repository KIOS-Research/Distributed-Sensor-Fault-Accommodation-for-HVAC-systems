% Test for the stability UNDER HEALTHY CONDITION ()

close all
clear all
clc

v=100; %1:1:50; % UNCERTAIN TERM
eo=50; % INITIAL TRACKING ERROR
A=-10; 
K=-9.5; %-100:100; % CONTROLLER GAIN
t=0:0.01:5; % TIME VARIABLE

f1=(A-K).*ones(length(t)); % A-K
f2=abs( v./(v+(eo-v).*exp((A-K).*t)) );
f=logical(f1<f2); % VALIDITY CHECK

figure
hold all
plot(t,f1,'color','blue','LineWidth',1.5)
plot(t,f2,':','color','red','LineWidth',1.5)
plot(t,f,'--','color','green','LineWidth',1.5)
h=legend('$A-K$','$|\frac{v}{v+(e0-v)exp((A-K)t)}|$','Valid')
    h.Interpreter='latex';
    h.FontSize=18;
    h.Orientation='vertical';
    h.Location='North';
ylim(gca,[min(min(f1(:)),min(min(f2(:)),min(f(:))))-1 max(max(f1(:)),max(max(f2(:)),max(f(:))))+1])