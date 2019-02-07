close all
clear all
clc

bias=1;
t=0:0.001:1;
xo=1;

a=4;
b=8*a;
c=55*bias;

F= @(t,x) -a*x - bias ;
x1=ode45(F,t,xo);

plot(x1.x,abs(x1.y))
s=x1.y;
hold on
F= @(t,x) -b*x + c ;
x1=ode45(F,t,xo);
plot(x1.x,abs(x1.y))
n= s<=x1.y;
plot(x1.x,logical(n))
ylim([-0.1 1.1])
% F= @(t,x) -1.5*x +bias 
% x1=ode45(F,0:0.001:10,xo)
% plot(x1.y,x1.x)