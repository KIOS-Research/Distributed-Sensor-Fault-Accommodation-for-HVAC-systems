t=0:0.1:100;

x=sin(t);
y=cos(2*t);

figure
hold on
plot(t,x)
plot(t,y)

plot3(t,x,y)

z=x.^2+y.^2
c=-2*(x.*y)

%isvalid(c<=z)

plot3(t,z,c)