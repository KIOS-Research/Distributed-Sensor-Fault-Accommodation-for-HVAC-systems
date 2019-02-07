function dx=system_f(t, x)

S=sat(-3.*x(2));
%dx=[-x(2); x(1)];

dx=[-x(2); 2.*x(1)+3.*x(2)+2.*S];
end