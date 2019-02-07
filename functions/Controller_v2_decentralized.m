function  U  = Controller_v2_decentralized(y, str_z, awz, h, p_air, Cp,...
    kst, Trefw, Ta, Cv, aw, Cw, Ustmax,...
    COPmax, DTmax, To, Tpl, Leng)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% global yy
% 
% yy=y;

%y(Leng+1)=


for k=1:(Leng)  
    Xref(k)= str_z(k).Tref;
    gain(k)= str_z(k).k; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g(k,k)=(((str_z(k).Umax * awz)) /(str_z(k).cz)) *...
        (Trefw - y(k));%(y(Leng+1) - y(k));
    itta(k)=(str_z(k).az/str_z(k).cz)*Ta...
        - ((h*str_z(k).Ad)/str_z(k).cz)*str_z(k).T1;
    sum4=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(k).connected))
        sum4=sum4+str_z(k).az_ij(c)/str_z(k).cz; %sum for A
        sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*...
            str_z(str_z(k).connected(c)).Tref;%str_z(k).Tref;
        for d=1:str_z(k).paths.TotalPaths
            sum3 = sum3 + sign( str_z(str_z(k).paths.ConnDoors(d)).Tref - y(k) )...
                * str_z(k).paths.Ad_ij(d)...
                * max( y(k) , str_z(str_z(k).paths.ConnDoors(d)).Tref ) *...
                sqrt(2 * (Cp - Cv) *...
                abs( str_z(str_z(k).paths.ConnDoors(d)).Tref - y(k)) );
        end
    end
    % str_z(k).Ad is the same this A_w_i
    A(k,k)=(h*str_z(k).Ad-(str_z(k).az))/str_z(k).cz - sum4;
    hh(k) = sum2 + ((p_air * Cp)/str_z(k).cz) * sum3;  
end
Uz=(-(y(1:Leng))*(A(1:Leng,1:Leng))-hh(1:Leng)-itta(1:Leng)...
    +gain.*(Xref-y(1:Leng)))/(g(1:Leng,1:Leng));

sum1=0;
for k=1:Leng
    sum1=sum1 + str_z(k).Umax * (y(Leng+1)-str_z(k).Tref)*Uz(k);
end

A(Leng+1,Leng+1)= -aw/Cw; %aw
g(Leng+1,Leng+1)= (Ustmax/Cw)* ...
    (1+(COPmax-1)*(1-((y(Leng+1)-To))/DTmax));
hh(Leng+1) = (awz/Cw)* sum1;%aw
itta(Leng+1) = (aw/Cw) * Tpl;%aw
%kst=10;
gain(Leng+1)=kst;
Xref(Leng+1)=Trefw;


U=(-y(1:Leng+1)*A-hh-itta+gain.*(Xref-y(1:Leng+1)))/g;

end

