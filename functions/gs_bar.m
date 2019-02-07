function bound=gs_bar(y,Ustmax,Cw,COPmax,DTmax,To,n_bar_s)

omega=Ustmax/Cw;
lambda=(COPmax-1)/DTmax;
kappa=DTmax+To;

if y <= kappa - n_bar_s
    bound=omega*lambda*n_bar_s;
elseif y > kappa + n_bar_s
    bound=0;
else
    bound= max(omega*lambda*abs(kappa -y -n_bar_s),omega*lambda*abs(kappa -y +n_bar_s));
end

end