%f=@(t,X)[-X(2); 2.*X(1)+3.*X(2)+2];%sat(-3.*X(2))
x1=linspace(-500,500,20);
x2=linspace(-500,500,20);
[x y]=meshgrid(x1,x2);
size(x)
size(y)
u=zeros(size(x));
v=zeros(size(x));
t=0;
for i=1:numel(x)
    Xprime=f(t,[x(i); y(i)]);
    u(i)=Xprime(1);
    v(i)=Xprime(2);
end
quiver(x,y,u,v,'k');figure(gcf)
xlabel('x_1')
ylabel('x_2')
axis tight equal;
hold on

f=@(t,X) system_f(t, X)


for i=-5:1:5
for x20=[-5 5]
    if (-5<i & i<5 & -5<x20 & x20<5)
        continue
    end
    [ts xs]=ode45(f,[0 1],[0;x20]);
    plot(xs(:,1),xs(:,2))
    plot(xs(1,1),xs(1,2),'ko')
    plot(xs(end,1),xs(end,2),'ks')
end
end
hold off