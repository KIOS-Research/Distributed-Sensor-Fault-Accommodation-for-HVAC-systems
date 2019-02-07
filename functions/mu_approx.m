clear all
close all
clc
%   x = -5:0.1:5;
%  y = -3.*sin(x);
% [x,y]=meshgrid(-2:0.1:2);
% t=0:0.1:10;
% 
% x=sin(t);
% y=cos(2*t);
   x = -1:0.1:1;
  y = -3.*(x);
  

figure
plot(x,sign(x),'LineWidth',2.5)
hold on
plot(x,tanh(3.*x),'LineWidth',2.5)
plot(x,tanh(10.*x),':','LineWidth',2.5)
set(gca,'FontSize',35)
h=legend('$sng(x)$','$\tanh(3x)$','$\tanh(10x)$');
    h.Interpreter='latex';
    h.FontSize=35;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',35) 

%%
for i=1:length(x)
    for j=1:length(y)
N(i,j)=sign(x(i)-y(j));
Ne(i,j)=tanh(3*(x(i)-y(j)));
    end
end


figure
% plot(x,N,'LineWidth',2.5)
surf(x,y,N,'FaceAlpha',0,'EdgeColor','b','FaceColor','b')
hold on
surf(x,y,Ne,'FaceAlpha',0,'EdgeColor','r','FaceColor','r')
% plot(x,Ne,':','LineWidth',2.5)
set(gca,'FontSize',18)
h=legend('$sgn(x-y)$','$\tanh(3(x-y))$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25) 
ylabel('$y$','Interpreter','latex','Fontsize',25)


%%

  
for i=1:length(x)
    for j=1:length(y)
M(i,j)=max(x(i),y(j));
    end
end


%M=max(x,y);
a=1;%a=1;
for i=1:length(x)
%    Me(i)=((x(i))*tanh(a*(x(i)-y(i)))) + ((y(i)+x(i))*tanh(a*(y(i)-x(i)))) - 2*((x(i))*tanh(a*(x(i)-y(i)))) ; % - ((y(i))*tanh(a*(x(i))));
    for j=1:length(y)
   Me(i,j)= ( x(i)*exp(a*x(i)) + y(j)*exp(a*y(j)) ) / ( exp(a*x(i)) + exp(a*y(j)) );
   Me(i,j)= log(exp(a*x(i)) + exp(a*y(j))); % max(x=0,y=0)=0 while log(exp(x=0) + exp(y=0)) is not equal to zero
    end   
% Me(i)=   ((x(i)-y(i))*tanh(a*(x(i)-y(i)))-(y(i))*tanh(a*(y(i)-x(i))))*tanh(a*(x(i)-y(i)))...
%     - ( (y(i)-x(i))*tanh(a*(y(i)-x(i))-y(i)*tanh(a*(x(i)-y(i)))) )*tanh(a*(x(i)-y(i))) ; % -(x(i))*tanh(x(i)-y(i));
end
%Me= ( x.*exp(a.*x) + y.*exp(a.*y) ) / ( exp(a.*x) + exp(a.*y) );

figure
%subplot(2,1,1)
% plot(x,M,'LineWidth',2.5)
 surf(x,y,M,'FaceAlpha',0,'EdgeColor','b','FaceColor','b')
hold on
surf(x,y,Me,'FaceAlpha',0,'EdgeColor','r','FaceColor','r')
 %plot(x,Me,':','LineWidth',2.5)
set(gca,'FontSize',18)
h=legend('$\max(x,y)$','$\frac{xe^{ax} + e^{ay}}{ e^{ax} + e^{ay}} $');%('M','Me');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25) 
ylabel('$y$','Interpreter','latex','Fontsize',25) 
% subplot(2,1,2)
% plot(y,abs(M(1,:)-Me(1,:)),'LineWidth',2.5)
% set(gca,'FontSize',18)


%%
%    x = -5:0.1:5;
%   y = -3.*sin(x);
%S=sqrt(abs(x-y));
ex=0.5;a=1;
% Se= (((x-y).^(0.5))./1).*tanh((x-y))+1;
for i=1:length(x)
    for j=1:length(y)
        S(i,j)=sqrt(abs(-x(i)+y(j)));
        Se(i,j)= ( ((x(i)*tanh(a*(x(i)-y(j)))-y(j)*tanh(a*(x(i)-y(j)))))^(0.5)  ).*tanh(a*(-x(i)+y(j)))...
            - ( ( x(i)*tanh(a*(x(i)-y(j))) - y(j)*tanh(a*(x(i)-y(j)))  )^(0.5) ).*tanh(a*(x(i)-y(j)));
%         ...
%             - (  ((y(j)*tanh(a*(y(j)-x(i)))+x(i)*tanh(a*(-y(i)+x(i)))))^(0.5) ).*tanh(a*(y(j)+x(i)));
        Sa(i,j)= ((x(i)-y(j))*exp(0.5*(x(i)-y(j)-1)) - (y(j)-x(i))*exp(0.5*(-x(i)+y(j))));
        Sa(i,j)=( ((x(i)*tanh(a*(x(i)-y(j)))-y(j)*tanh(a*(x(i)-y(j)))))^(0.5)  );
    end
end
%         S=sqrt(abs(x-y));
%         Se= ( ((x.*tanh(a.*(x-y))-y.*tanh(a.*(x-y)))).^(0.5)  ).*tanh(a.*(x-y)) ...
%            + (  ((y.*tanh(a.*(y-x))-x.*tanh(a.*(y-x)))).^(0.5) );
figure
surf(x,y,S,'FaceAlpha',0,'EdgeColor','b','FaceColor','b')
hold on
% surf(x,y,Se,'FaceAlpha',0,'EdgeColor','r')
surf(x,y,Sa,'FaceAlpha',0,'EdgeColor','r','FaceColor','r')
% plot(x,S,'LineWidth',2.5)
% hold on
% plot(x,Se,':','LineWidth',2.5)
set(gca,'FontSize',18)
h=legend('$ {\sqrt{ \left| x-y \right| }}$','$ {\sqrt{  (x-y)tanh(a(x-y)) }}$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25) 
%%

figure
surf(x,y,N.*M.*S,'FaceAlpha',0,'EdgeColor','b','FaceColor','b')
hold on
surf(x,y,Ne.*Me.*Sa,'FaceAlpha',0,'EdgeColor','r','FaceColor','r')
% plot(x,N.*M.*S,'LineWidth',2.5)
% hold on
% plot(x,Ne.*Me.*Se,':','LineWidth',2.5)
set(gca,'FontSize',18)
h=legend('$\mu(.)$','$\hat{\mu}(.)$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25) 



%%
% x=-2:0.1:2;
% y=-2:0.1:2;
% for i=1:length(x)
%     for j=1:length(y)
% z(i,j)=(x(i).*exp(x(i)) + y(j).*exp(y(j))) ./ (exp(x(i))+exp(y(j)));
% zm(i,j)=max(x(i),y(j));
%     end
% end

[X,Y]=meshgrid(-2:0.2:2);
z=(X.*exp(X) + Y.*exp(Y)) ./ (exp(X)+exp(Y));
zm=max(X,Y);

% plot3(x,y,z)
surf(X,Y,z,'FaceAlpha',0.6,'EdgeColor','b');
hold on
surf(X,Y,zm,'FaceAlpha',0.6,'EdgeColor','r');


%%
z1=sign(x);
for i=1:length(x)
    k1(i)=max(x(i));
end
m1=sqrt(abs(x));

z2=tanh(3.*x);
for i=1:length(x)
    k2(i)=x(i);%*tanh(x(i));%-x(i)*tanh(x(i));
end
m2= (x./4).*tanh(x)+1;
% for i=1:length(x)
%     for j=1:length(y)
% z(i,j)=tanh(x(i)-y(j));
% k(i,j)=x(i)*tanh(x(i)-y(j))+y(j)*tanh(y(j)-x(i));
%     end
% end

figure
plot(x,z1,'LineWidth',2.5)
hold on
plot(x,z2,'LineWidth',2.5,'LineStyle',':')
set(gca,'FontSize',18)
h=legend('$sgn(x)$','$\tanh(3x)$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25) 

figure
plot(x,k1,'LineWidth',2.5)
hold on
plot(x,k2,'LineWidth',2.5,'LineStyle',':')
set(gca,'FontSize',18)
h=legend('$\max(x)$','$x$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25)  

figure
plot(x,m1,'LineWidth',2.5)
hold on
plot(x,m2,'LineWidth',2.5,'LineStyle',':')
set(gca,'FontSize',18)
h=legend('$ {\sqrt{ \left| x \right| }}$','$\left(\frac{x}{4}\tanh(x)+1\right)$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25)  

% compute polynomial
g=0.00*x.^(3)+0.34*x.^(2)+0.001*x;

figure
plot(x,z1.*k1.*m1,'LineWidth',2.5)
hold on
plot(x,z2.*k2.*m2,'LineWidth',2.5,'LineStyle',':')
hold on
plot(x,g,'LineWidth',2.5,'LineStyle','--')
set(gca,'FontSize',18)
h=legend('$sgn(x)\max(x) {\sqrt{ \left| x \right| }}$','$\tanh(3x)x\left(\frac{x}{4}\tanh(x)+1\right)$','$ax^2+bx+c$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x$','Interpreter','latex','Fontsize',25)    
% figure
% plot3(x,y,z)
% figure
% plot3(x,y,k)

%% plot derivative

syms x1 x2

g= tanh(x2-x1)*log(exp(x1)+exp(x2))*sqrt((x2-x1)*tanh(x2-x1));

f=jacobian(g,[x1 x2]);

[x1,x2]=meshgrid(-40:2:40);


for i=1:length(x1)
    for j=1:length(x2)
       G(i,j)=tanh(x2(j)-x1(i))*log(exp(x1(i))+exp(x2(j)))*sqrt((x2(j)-x1(i))*tanh(x2(j)-x1(i))) ;
       F1(i,j)=log(exp(x1(i)) + exp(x2(j)))*(tanh(x1(i) - x2(j))*(x1(i) - x2(j)))^(1/2)*(tanh(x1(i) - x2(j))^2 - 1)...
           - (tanh(x1(i) - x2(j))*log(exp(x1(i)) + exp(x2(j)))*(tanh(x1(i) - x2(j)) - (x1(i) - x2(j))*(tanh(x1(i) - x2(j))^2 - 1)))...
           /(2*(tanh(x1(i) - x2(j))*(x1(i) - x2(j)))^(1/2)) - (tanh(x1(i) - x2(j))*exp(x1(i))*(tanh(x1(i) - x2(j))*(x1(i) - x2(j)))^(1/2))/(exp(x1(i)) + exp(x2(j)));
    end
end


figure
%surf(x1,x2,G,'FaceAlpha',0,'EdgeColor','b','FaceColor','b')
hold on
surf(x1,x2,-(str_z(1).az_ij(1)/str_z(1).cz)/((sqrt(2 * (Cp - Cv))*str_z(1).paths.Ad_ij)/str_z(1).cz)...
    *ones(length(x1),length(x2)),'FaceAlpha',0,'EdgeColor','g','FaceColor','g')
hold on
surf(x1,x2,F1,'FaceAlpha',0,'EdgeColor','r','FaceColor','r')
set(gca,'FontSize',18)
h=legend('$-a_{12}$','$ \phi_{12}$');
    h.Interpreter='latex';
    h.FontSize=25;
    h.Orientation='vertical';
    h.Location='North';
xlabel('$x_1$','Interpreter','latex','Fontsize',25) 
ylabel('$x_2$','Interpreter','latex','Fontsize',25) 


