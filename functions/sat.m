function Y=sat(s)
% sat is the saturation function with unit limits and unit slope.
if s>1
Y=1;
elseif s<-1 
Y=-1;
else 
Y=s;
end