function  U  = Controller_v5(y, x, str_z, awz, h, p_air, Cp,...
    kst, Trefw, Ta, Cv, aw, Cw, Ustmax,...
    COPmax, DTmax, To, Tpl, AC, Leng)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% global yy
% 
% yy=y;

% global DD KK yv AC



for k=1:(Leng)
    
    Xref(k)= str_z(k).Tref;
    gain(k)= str_z(k).k; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    
%     if DD(k)==1 && KK(k)==0 % (LOCAL AGENT DETECTED A FAULT) e.g., k=2
%         y(k)=str_z(k).Tref;% y(k)+(x(k+Leng+1)-y(k)); 
%     elseif KK(k)==1
%         y(k)=yv(k);
%     end
    if AC(k)==1
        y(k)=y(k)-x(3*Leng+3+k);
    end
    
    g(k,k)=(((str_z(k).Umax * awz)) /(str_z(k).cz)) *...
        (y(Leng+1) - y(k));
    itta(k)=(str_z(k).az/str_z(k).cz)*Ta...
        - ((h*str_z(k).Ad)/str_z(k).cz)*str_z(k).T1;
    sum4=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(k).connected))
        sum4=sum4+str_z(k).az_ij(c)/str_z(k).cz; %sum for A
        sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*...
            y(str_z(k).connected(c));

%         if  DD(k)==1 % (LOCAL AGENT DETECTED A FAULT) e.g., k=2
%             sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*...
%                       x(str_z(k).connected(c)+Leng+1); %str_z(str_z(k).connected(c)).Tref; %x(str_z(k).connected(c)+Leng+1); %y(str_z(k).connected(c));
% %               sum2=sum2+(str_z(k).az_ij(c)/str_z(k).cz)*...
% %                       x(str_z(k).connected(c)+Leng+1); %str_z(str_z(k).connected(c)).Tref;%
%               
%         end
        
        for d=1:str_z(k).paths.TotalPaths
            sum3 = sum3 + sign( y(str_z(k).paths.ConnDoors(d))-y(k))...
                * str_z(k).paths.Ad_ij(d)...
                * max(y(k),y(str_z(k).paths.ConnDoors(d))) *...
                sqrt(2 * (Cp - Cv) *...
                abs(y(str_z(k).paths.ConnDoors(d))-y(k)));
%             if DD(k)==1 %|| DD(str_z(k).paths.ConnDoors(d))==1 % (NEIGHBOR AGENT DETECTED A FAULT) e.g., k=1,3 
%                   sum3 =sum3 + sign( y(str_z(k).paths.ConnDoors(d)) - x(k+Leng+1) )...
%                    * str_z(k).paths.Ad_ij(d)...
%                    * max( x(k+Leng+1) , y(str_z(k).paths.ConnDoors(d)) ) *...
%                    sqrt(2 * (Cp - Cv) *...
%                    abs( y(str_z(k).paths.ConnDoors(d)) - x(k+Leng+1) ) );
%             end
%             if DD(k)==1 % (LOCAL AGENT DETECTED A FAULT) e.g., k=2
% %                  sum3 =sum3 + sign( y(str_z(k).paths.ConnDoors(d)) - str_z(k).Tref )...
% %                    * str_z(k).paths.Ad_ij(d)...
% %                    * max( str_z(k).Tref , y(str_z(k).paths.ConnDoors(d)) ) *...
% %                    sqrt(2 * (Cp - Cv) *...
% %                    abs( y(str_z(k).paths.ConnDoors(d)) - str_z(k).Tref ) );
%                  sum3 = sum3 + sign( x(str_z(k).paths.ConnDoors(d)+Leng+1) - y(k) )...
%                     * str_z(k).paths.Ad_ij(d)...
%                     * max( y(k) , x(str_z(k).paths.ConnDoors(d)+Leng+1) ) *...
%                     sqrt(2 * (Cp - Cv) *...
%                     abs( x(str_z(k).paths.ConnDoors(d)+Leng+1) - y(k)) );
%             end
        end
    end
    % str_z(k).Ad is the same this A_w_i
    A(k,k)=(h*str_z(k).Ad-(str_z(k).az))/str_z(k).cz - sum4;
    hh(k) = sum2 + ((p_air * Cp)/str_z(k).cz) * sum3;
    
end

% Uz=(-(y(1:Leng))*(A(1:Leng,1:Leng))-hh(1:Leng)-itta(1:Leng)...
%     +gain.*(Xref-y(1:Leng)))/(g(1:Leng,1:Leng));
Uz=(-Xref*(A(1:Leng,1:Leng))-hh(1:Leng)-itta(1:Leng)...
    -gain.*(y(1:Leng)-Xref))/(g(1:Leng,1:Leng));

sum1=0;
for k=1:Leng
    sum1=sum1 + str_z(k).Umax * (y(Leng+1)-y(k))*Uz(k);
end

A(Leng+1,Leng+1)= -aw/Cw; %aw
g(Leng+1,Leng+1)= (Ustmax/Cw)* ...
    (1+(COPmax-1)*(1-((y(Leng+1)-To))/DTmax));
hh(Leng+1) = (awz/Cw)* sum1;%aw
itta(Leng+1) = (aw/Cw) * Tpl;%aw
%kst=10;
gain(Leng+1)=kst;
Xref(Leng+1)=Trefw;


% U=(-y(1:Leng+1)*A-hh-itta+gain.*(Xref-y(1:Leng+1)))/g;
U=(-Xref*A-hh-itta-gain.*(y(1:Leng+1)-Xref))/g;
end

