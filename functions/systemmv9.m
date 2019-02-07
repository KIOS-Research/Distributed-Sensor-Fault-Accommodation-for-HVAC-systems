
function dx=systemmv9(t, x, time, step, tspan, tt, nn, str_z, k, con,...
    Az ,Umax, cz ,az, Ad, T1, h, awz, Cw, Ustmax,...
    COPmax, To, DTmax, Cp,...
    Cv, p_air, Ta, Tpl, aw, maxst,...
    minst, kst, Trefw, L_s, n_bar_s,...
    p_s, lambda_s, x_bar_s, r_bar_s,...
    K, Leng, amp_out, amp_s, f_out, f_s, F_time_s, F_value_s,DET,Df)

n=zeros(1,Leng+1);
y=zeros(1,Leng+1);
yV=zeros(1,Leng+1);
f=zeros(1,Leng+1);
fa=zeros(1,Leng+1);
DD=zeros(1,Leng+1);
TD=zeros(1,Leng+1);
AC=zeros(1,Leng+1);
DDa=zeros(1,Leng);
DDs=zeros(1,Leng);
gu=zeros(1,Leng+1);
gubar=zeros(1,Leng);
Dg=0;
gbar=0;


global TD TDa TDs TLs H DD U U_time PP KK W yv AC%sum88
%sumt mm
% Sensor Fault

%global mm sumt n

tt=find(tspan==t); % the number of solution 


if t==0
    y(1:Leng+1)=x(1:Leng+1);  %x(3*Leng+4:4*Leng+4);
else
    
    for k=1:Leng
        % injecting the sensor faults
        if t<str_z(k).F_time
            f(k)=0;
        else
            f(k)=str_z(k).F_value*(1-exp(-K*t));
        end
        % injecting the actuator faults
        if t<str_z(k).Fa_time
            fa(k)=0;
        else
            fa(k)=str_z(k).Fa_value*(1-exp(-K*t));
        end
%         % reducing the samples of noise
%         if tt ==  nn 
%             nn=2*nn;
%             n1=str_z(k).max+2*str_z(k).min*rand(1); %%%%%% noise at sensor of zones
%         else
%             n1=n(k);   
%         end
%         n(k)=n1; %str_z(k).max+2*str_z(k).min*rand(1); %n1;
        n(k)=str_z(k).max+2*str_z(k).min*rand(1);
        %n(k)=filter(2,[1 2],n(k));
        y(k)=x(k)+n(k)+f(k);
    end

end


%noise for the storage tank

%         % reducing the samples of noise
%         if tt ==  nn 
%             nn=2*nn;
%             n(Leng+1)= maxst+2*minst*rand(1); %%%%%% noise at sensor of zones
%         else
%             n(Leng+1)=n(Leng+1);   
%         end

        if t<F_time_s
            f(Leng+1)=0;
        else
            f(Leng+1)=F_value_s*(1-exp(-K*t));
        end
n(Leng+1)= maxst+2*minst*rand(1);
y(Leng+1) = x(Leng+1) + f(Leng+1) + n(Leng+1);
%fa(Leng+1)=0;
%global mm U tt D E Ebar Ea yy
%tt{mm}=t;
% mm
% tt{mm}
%%%%%%%%%%%%%%%%%
% if t==0
%         T=[];
% end
% T(end+1)=t;

% SW=0;
% fval=[];
% for k=1:Leng
% %  if DD(k) == 1
% %   str_z(k).k =K_f(k);
% %  end
%     if TD(k)<= t
%         SW=1;
% %         break;
%     end
%     k_prev=str_z(k).k;
%     if t==0
%         U(k)=0;
%         U_time=[];
%         T=[];
%         sum88(1,k)=0;
%         flag=0;
%         fval(k)=2;
% 
%     else
%      [str_z(k).k,fval(k)]=Optimal_K_u(Leng,k,t,y,str_z,h,awz,r_bar_s, n_bar_s,...
%         Cp, Cv,p_air,Ta,U_time(:,k),T,k_prev);   
%     end
% end

%%%%%%%%%%%%%%
str_z(2).k
t
% if SW==0
U = Controller_v5(y, x, str_z, awz, h, p_air,...
    Cp, kst, Trefw, Ta, Cv, aw, Cw, Ustmax,...
    COPmax, DTmax, To, Tpl, AC, Leng);
U_time(end+1,:)=U(1:Leng);

for i=1:(length(U)-1) %(length(U)-1)
if U(i) >= 1
    U(i)=1;
elseif U(i) <= 0
    U(i)=0;
end
end
if U(length(U)) <= -1
    U(length(U))=-1;
elseif U(length(U)) >= 0
    U(length(U))=0;
end
% else
% U = Controller_v4_decentralized(y, x, str_z, awz, h, p_air,...
%     Cp, kst, Trefw, Ta, Cv, aw, Cw, Ustmax,...
%     COPmax, DTmax, To, Tpl, Leng);
%end

sum4=0;
sum2=0;
sum3=0;



%%%%%%%%%%%%%% With matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     g(k)=((Umax(k) * awz) ./ cz(k)) .* (x(Leng+1) - x(k))';
%     L(k)=str_z(k).L;
%
%     sum1= Umax(k) * (x(Leng+1)*eye(Leng) -diag(x(k))) * U{mm}(k)';
%         itta(k)=(az(k) ./ cz(k))*Ta - ...
%             ((h*Ad(k))./ cz(k)) .* T1(k);
%%%%%%%%%%%%%% With matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%


sum1=0;
for i=1:(Leng)
    %i
    g(i,i)=(((str_z(i).Umax * awz)) /(str_z(i).cz)) *...
        (x(Leng+1) - x(i));
    L(i,i)=str_z(i).L;
    sum1=sum1 + str_z(i).Umax * (x(Leng+1)-x(i))*(U(i));
    itta(i)=(str_z(i).az/str_z(i).cz)*Ta...
        - ((h*str_z(i).Ad)/str_z(i).cz)*str_z(i).T1;
    sum4=0;
    sum2=0;
    sum3=0;

    for c=1:(length(str_z(i).connected))     %sum for A
        sum4=sum4+str_z(i).az_ij(c)/str_z(i).cz;
        sum2=sum2+(str_z(i).az_ij(c)/str_z(i).cz)*...
            x(str_z(i).connected(c));
        
        for d=1:str_z(i).paths.TotalPaths
            sum3 = sum3 + sign( x(str_z(i).paths.ConnDoors(d)) - x(i) )*...
                str_z(i).paths.Ad_ij(d) *...
                max(x(i),x(str_z(i).paths.ConnDoors(d))) *...
                sqrt(2 * (Cp - Cv) *...
                abs(x(str_z(i).paths.ConnDoors(d)) - x(i)));
            % str_z(i).paths.Ad_ij(d) * x(i) *...
            
        end
    end
    A(i,i)=(h*str_z(i).Ad-(str_z(i).az)) /str_z(i).cz - sum4;
    hh(i) = sum2 + ((p_air * Cp) / str_z(i).cz) * sum3;
    rho(i) =(str_z(i).az/str_z(i).cz)*amp_out*sin(f_out*t); 

end
% C4=[1 0 0 0];
% A4=[A(1,1) str_z(1).az_ij(1)/str_z(1).cz 0 0;
%    str_z(1).az_ij(1)/str_z(2).cz A(2,2) str_z(2).az_ij(2)/str_z(2).cz 0;
%    0 str_z(2).az_ij(2)/str_z(3).cz A(3,3) str_z(3).az_ij(1)/str_z(3).cz;
%    0 0 str_z(3).az_ij(1)/str_z(4).cz A(4,4)]
% Q4=[C4;C4*A4;C4*A4*A4;C4*A4*A4*A4]
% 
% C3=[1 0 0]
% A3=[A(1,1) str_z(1).az_ij(1)/str_z(1).cz 0;
%    str_z(1).az_ij(1)/str_z(2).cz A(2,2) str_z(2).az_ij(2)/str_z(2).cz;
%    0 str_z(2).az_ij(2)/str_z(3).cz A(3,3)]
% Q3=[C3;C3*A3;C3*A3*A3]
% inv(Q3)


A(Leng+1,Leng+1)= -aw/Cw;%aw
if x(Leng+1)-To <= DTmax
    Ps=1+ (COPmax-1)*(1-(x(Leng+1)-To)/DTmax);
else
    Ps=1;
end
g(Leng+1,Leng+1)= (Ustmax/Cw)*Ps;
%...    (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));

glim=(Ustmax/Cw);

Dg=abs(g(Leng+1,Leng+1)-glim);

gbar=(Ustmax*(COPmax-1))/(Cw*DTmax);

hh(Leng+1) = (awz/Cw)* sum1; %aw
itta(Leng+1) = (aw/Cw) * Tpl;%aw
rho(Leng+1) = (awz/Cw)*amp_s*sin(f_s*t);

%system equations (zones & storage tank)
dx1=x(1:Leng+1)'*A(1:Leng+1,1:Leng+1)...
    + (U(1:Leng+1)+ fa)*g(1:Leng+1,1:Leng+1) + hh(1:Leng+1)...
    + itta(1:Leng+1) + rho(1:Leng+1);

gu1=(U(1:Leng+1)+ fa)*g(1:Leng+1,1:Leng+1);

clear i m c



%%%%%

% Detection observer
sum1=0;
for i=Leng+2:2*Leng+1
 
    m=i-(Leng+1);
%     yy(m)=y(m);
%     
%     if DD(m)==1 && KK(m)==0 % (LOCAL AGENT DETECTED A FAULT) e.g., k=2
%         yy(m)=str_z(m).Tref;% y(k)+(x(k+Leng+1)-y(k)); 
%     elseif KK(m)==1
%         yy(m)=yv(m);
%     end
%     
    g(i,i)=(((str_z(m).Umax * awz)) /(str_z(m).cz)) *...
        (y(Leng+1) - y(m));
%     if DD(m)==1
%         g(i,i)=(((str_z(m).Umax * awz)) /(str_z(m).cz)) *...
%         (y(Leng+1) - str_z(m).Tref);
%     end
    
    sum1=sum1 + str_z(m).Umax * (x(Leng+1)-x(i))*U(m);
    itta(i)=(str_z(m).az/str_z(m).cz)*Ta...
        - ((h*str_z(m).Ad)/str_z(m).cz)*str_z(m).T1;
    sum4=0;
    sum2=0;
    sum3=0;
    for c=1:(length(str_z(m).connected))     %sum for A
        sum4=sum4+str_z(m).az_ij(c)/str_z(m).cz;
        sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)...
            *y(str_z(m).connected(c));
%         if  DD(str_z(m).connected(c))==1 %(NEIGHBOR AGENT DETECTED A FAULT) e.g., k=1,3 
%             sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)*...
%                    str_z(str_z(m).connected(c)).Tref; % x(str_z(m).connected(c)+Leng+1); %
% %               sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)*...
% %                       x(str_z(m).connected(c)+Leng+1); %str_z(str_z(m).connected(c)).Tref;%
%         end
%         if  DD(m)==1 % (LOCAL AGENT DETECTED A FAULT) e.g., m=2
%             sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)*...
%                        str_z(str_z(m).connected(c)).Tref; %x(str_z(m).connected(c)+Leng+1); %str_z(str_z(m).connected(c)).Tref; %y(str_z(m).connected(c)); %
% %               sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)*...
% %                       x(str_z(m).connected(c)+Leng+1); %str_z(str_z(m).connected(c)).Tref;%
%         end
        
        for d=1:str_z(m).paths.TotalPaths
            sum3 = sum3 +...
                sign( y(str_z(m).paths.ConnDoors(d)) - y(m) ) *...
                str_z(m).paths.Ad_ij(d) *...
                max(y(m),y(str_z(m).paths.ConnDoors(d))) *...
                sqrt(2 * (Cp - Cv) * ...
                abs(y(str_z(m).paths.ConnDoors(d)) - y(m)));
            % str_z(m).paths.Ad_ij(d) * x(i) *...
%             if DD(str_z(m).paths.ConnDoors(d))==1 % (NEIGHBOR AGENT DETECTED A FAULT) e.g., m=1,3 
% %                   sum3 =sum3 + sign( x(str_z(m).paths.ConnDoors(d)+Leng+1) - x(m+Leng+1) )...
% %                    * str_z(m).paths.Ad_ij(d)...
% %                    * max( x(m+Leng+1) , x(str_z(m).paths.ConnDoors(d)+Leng+1) ) *...
% %                    sqrt(2 * (Cp - Cv) *...
% %                    abs( x(str_z(m).paths.ConnDoors(d)+Leng+1) - x(m+Leng+1) ) );
%                   sum3 =sum3 + sign( str_z(str_z(m).paths.ConnDoors(d)).Tref - str_z(m).Tref )...
%                    * str_z(m).paths.Ad_ij(d)...
%                    * max( str_z(m).Tref , str_z(str_z(m).paths.ConnDoors(d)).Tref ) *...
%                    sqrt(2 * (Cp - Cv) *...
%                    abs( str_z(str_z(m).paths.ConnDoors(d)).Tref - str_z(m).Tref ) );
%             end
%             if DD(m)==1 % (LOCAL AGENT DETECTED A FAULT) e.g., m=2
% %                  sum3 =sum3 + sign( y(str_z(m).paths.ConnDoors(d)) - str_z(m).Tref )...
% %                    * str_z(m).paths.Ad_ij(d)...
% %                    * max( str_z(m).Tref , y(str_z(m).paths.ConnDoors(d)) ) *...
% %                    sqrt(2 * (Cp - Cv) *...
% %                    abs( y(str_z(m).paths.ConnDoors(d)) - str_z(m).Tref ) );
%                   sum3 =sum3 + sign( str_z(str_z(m).paths.ConnDoors(d)).Tref - str_z(m).Tref )...
%                    * str_z(m).paths.Ad_ij(d)...
%                    * max( str_z(m).Tref , str_z(str_z(m).paths.ConnDoors(d)).Tref ) *...
%                    sqrt(2 * (Cp - Cv) *...
%                    abs( str_z(str_z(m).paths.ConnDoors(d)).Tref - str_z(m).Tref ) );
%                    
%                y(m)=str_z(m).Tref;
%                
%             end
        end
    end
    A(i,i)=(h*str_z(m).Ad-(str_z(m).az))/str_z(m).cz - sum4;
    hh(i) = sum2 + ((p_air * Cp)/str_z(m).cz) * sum3;
end

% Tref(1:Leng)=str_z(1:Leng).Tref;
% Tref(Leng+1)=Trefw;

A(2*Leng+2,2*Leng+2)= -aw/Cw;%aw
g(2*Leng+2,2*Leng+2)= (Ustmax/Cw)*...
    (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));
hh(2*Leng+2) = (awz/Cw)* sum1; %aw
itta(2*Leng+2) = (aw/Cw) * Tpl;%aw

% storage tank observer gain
L(Leng+1,Leng+1)=L_s;

% Estimator (zones & storage tank)

LEst=Leng+2:2*Leng+2;
dxe= x(LEst)'*A(LEst,LEst)...
    + U*g(LEst,LEst)...
    + hh(LEst) ...
    + itta(LEst)...
    + (y(1:Leng+1)-x(LEst)')*L;


gu2=U*g(LEst,LEst);

gu=gu1-gu2;

% Detection residual
E(1:Leng+1)=y(1:Leng+1)-(x(LEst))';

%%
% Detection threshold
sum5=0;
for l=1:Leng
    sum5=sum5 + str_z(l).Umax * (n_bar_s + str_z(l).n_bar)...
        * abs(U(l));
    sum6=0;
    
    for d=1:(length(str_z(l).connected))     %sum for A
        
        sum6 = sum6 + str_z(l).az_ij(d)...
            * (str_z(l).n_bar  + str_z(d).n_bar);
        sum7=0;
        for p=1:str_z(l).paths.TotalPaths
            %Dmu=bound(y(l),y(str_z(l).paths.ConnDoors(b)),...
            %    str_z(l).n_bar,str_z(str_z(l).paths.ConnDoors(b)).n_bar);
         c=sign(y(str_z(l).paths.ConnDoors(p))-y(l))...
            .*max(y(str_z(l).paths.ConnDoors(p)),y(l))...
            .*sqrt(abs(y(str_z(l).paths.ConnDoors(p))-y(l)));


    if y(str_z(l).paths.ConnDoors(p)) -y(l) >...
            str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        a2 =y(str_z(l).paths.ConnDoors(p))...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        b2 =y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        
    elseif y(str_z(l).paths.ConnDoors(p)) -y(l) <...
            -str_z(l).n_bar-str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar)); 
        a2 =-y(l) -str_z(l).n_bar;
        b2 =-y(l) +str_z(l).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    else % |y(str_z(l).paths.ConnDoors(b))-y(l)|<str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(b)).n_bar
     
        a1 =0;
        b1 =sqrt(max(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar),abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar)));
        a2 =min([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        b2 =max([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    end

        Dmu =max(abs(a -c ),abs(b -c ));
            
        sum7 = sum7 + str_z(l).paths.Ad_ij(p) * Dmu;  % h_bar(l);
        
        end
        
        
        e1= abs( str_z(l).L ) * str_z(l).n_bar + str_z(l).r_bar...
            + ((str_z(l).Umax * awz) / str_z(l).cz) * ...
            (str_z(l).n_bar + n_bar_s) * abs(U(l))...
            + (1/str_z(l).cz) * sum6...
            + ((p_air * Cp)/str_z(l).cz) *...
            sqrt( 2 * abs(Cp - Cv)) * sum7;
        
        gubar(l)=((str_z(l).Umax * awz) / str_z(l).cz) * ...
            (str_z(l).n_bar + n_bar_s) * abs(U(l));
    end
    
    funi=@(t) str_z(l).p*exp(-str_z(l).lambda * t)*e1;
    
    %     sum1=0;
    %     for i=1:tt
    %         sum1=sum1 + str_z(l).lambda^(tt-i)*e1{i};
    %     end
    
    % threshold of zone i
    Ebar(l) = str_z(l).p * exp(-str_z(l).lambda * t)...
        *str_z(l).x_bar...
        +str_z(l).n_bar +integral(funi,0,t);
    % + str_z(l).p*sum1;
    % +integral(funi,0,tt);
    %              +lsim(TF1,e1,t);
    %,'RelTol',1e-4,'AbsTol',1e-4);
end


%%
%TF =tf([0 p_s] , [1  lambda_s]);

% es= abs((Ustmax * (COPmax-1))/(Cw * DTmax)) * n_bar_s *...
%     abs( U(Leng+1) ) + abs(L_s) * n_bar_s;

es= gs_bar(y(Leng+1),Ustmax,Cw,COPmax,DTmax,To,n_bar_s)*...
    abs( U(Leng+1) ) + abs(L_s) * n_bar_s;

funs=@(t) p_s*exp(-lambda_s * t)*es;

% sum1=0;
% for i=1:tt
%     sum1=sum1 + lambda_s^(tt-i)*es{i};
% end

% threshold of storage tank
Ebar(Leng+1) = p_s * exp(-lambda_s * t) * x_bar_s ...
    +(awz/Cw) *sum5 + r_bar_s...
    +n_bar_s... 
    +integral(funs,0,t);
%+ (aw/(Cw*1000)) *sum5 + r_bar_s...
% ,'RelTol',1e-4,'AbsTol',1e-4) ...
%+integral(funs,0,tt)...
%+ p_s*sum1..


% Decision

for k=1:Leng+1
    if abs(E(k))>Ebar(k)
        DD(k)=1;
        TD(k)=min(TD(k),t);
        
%         PP=PP+1;
%         if PP<10   
%             H(PP,l)=abs(Es(l));
%         else
%             PP=0;
%         end
% 
%         if abs(Es(l)) <= Ebars(l)
%             TLs(l)=min(TD(l),t);
%         end
        
    else
        DD(k)=0;
    end
    
    if t>= TD(k)
        DD(k)=1;
    end
        
end


%%%%%%%%%%%%% without adaptive estimation

I=2*Leng+3:3*Leng+3;

 
    
sum1=0;
for i=2*Leng+3:3*Leng+2 %Leng+2:2*Leng+1
    
    m=i-(2*Leng+2); %i-(Leng+1);
%     if t==TD(m)
%      x(I(m))=y(m);
%     end
    
    if t>=TD(m)+0.2 %0.4
%           
        AC(m)=1;

        g(i,i)=(((str_z(m).Umax * awz)) /(str_z(m).cz)) *...
            (y(Leng+1) - x(i));
        G(m)=(((str_z(m).Umax * awz)) /(str_z(m).cz));
        sum1=sum1 + str_z(m).Umax * (y(Leng+1)-x(i))*U(m);
        itta(i)=(str_z(m).az/str_z(m).cz)*Ta...
            - ((h*str_z(m).Ad)/str_z(m).cz)*str_z(m).T1;
        sum4=0;
        sum2=0;
        sum3=0;
        for c=1:(length(str_z(m).connected))     %sum for A
            sum4=sum4+str_z(m).az_ij(c)/str_z(m).cz;
            sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)...
                *y(str_z(m).connected(c)); %x(I(str_z(m).connected(c)));

            for d=1:str_z(m).paths.TotalPaths
                sum3 = sum3 +...
                    sign( y(str_z(m).paths.ConnDoors(d)) - x(i) ) *...
                    str_z(m).paths.Ad_ij(d) *...
                    max(x(i),y(str_z(m).paths.ConnDoors(d)) ) *...
                    sqrt(2 * (Cp - Cv) * ...
                    abs(y(str_z(m).paths.ConnDoors(d)) - x(i)) );
               %x(I(str_z(m).paths.ConnDoors(d)))
            end
        end
        A(i,i)=(h*str_z(m).Ad-(str_z(m).az))/str_z(m).cz - sum4;
        hh(i) = sum2 + ((p_air * Cp)/str_z(m).cz) * sum3;

%      else
%        dxes(m)=zeros(1,Leng+1);
%        Es=zeros(1,Leng+1);
     end
    
end
AC(Leng+1)=0;

A(3*Leng+3,3*Leng+3)= -aw/Cw;%aw
g(3*Leng+3,3*Leng+3)= (Ustmax/Cw)*...
    (1+(COPmax-1)*(1-((y(Leng+1)-To))/DTmax));
G(Leng+1)=0;
hh(3*Leng+3) = (awz/Cw)* sum1; %aw
itta(3*Leng+3) = (aw/Cw) * Tpl;%aw

L_Est=2*Leng+3:3*Leng+3;
dxes= (x(L_Est)'*A(L_Est,L_Est)...
    + U*g(L_Est,L_Est)...
    + hh(L_Est) ...
    + itta(L_Est)).*AC ; %...
   % + (y(1:Leng+1)-x(LEst)')*L;

% Detection residual
Es(1:Leng+1)=(y(1:Leng+1)-(x(L_Est))').*AC;

% moving average
% windowSize=5;
% b=(1/windowSize)*ones(1,windowSize);
% a=1;

% Lowpass Butterworth Transfer Function
% fc = 300;
% fs = 1000;
% [b,a] = butter(1,fc/(fs/2));

% b=[0 1];
% a=[1 1];
% fhat=filter(b,a,Es);

% d = fdesign.lowpass('Fp,Fst,Ap,Ast',3,5,0.5,40,100);
%    Hd = design(d,'equiripple');
%    
%   Hd = tf(4 , [1 4]); % -5 /(s-1)
%   sys=ss(Hd)
fhat=Es; %filter(2, [1 2],Es);

Lf=3*Leng+4:4*(Leng)+4;
%dfhats=( -30.*x(Lf)' + 30.*Es).*AC;

kk(1:Leng)=str_z(1:Leng).k;
kk(Leng+1)=kst;

for i=1:Leng
R(i)=str_z(i).Tref;
end
R(Leng+1)=Trefw;

%dfhats=2*(kk - G.*U).*(y-R).*AC;
dfhats= 4*(G.*U-kk).*(y-R-x(Lf)').*AC;

%fhat=filter(5,[5 1],Es,1); (5/5+z-1) 
z=y-x(Lf)';
s=y-R-(x(Lf)');


if DET==0

%%
% healthy tracking error threshold
% sum5=0;
for l=1:Leng
%     sum5=sum5 + str_z(l).Umax * (n_bar_s + str_z(l).n_bar)...
%         * abs(U(l));
    sum6=0;
    
    for d=1:(length(str_z(l).connected))     %sum for A
        
        sum6 = sum6 + str_z(l).az_ij(d)...
            * (str_z(d).n_bar);
        sum7=0;
        sum8=0;
        for p=1:str_z(l).paths.TotalPaths
            %Dmu=bound(y(l),y(str_z(l).paths.ConnDoors(b)),...
            %    str_z(l).n_bar,str_z(str_z(l).paths.ConnDoors(b)).n_bar);
         c=sign(y(str_z(l).paths.ConnDoors(p))-y(l))...
            .*max(y(str_z(l).paths.ConnDoors(p)),y(l))...
            .*sqrt(abs(y(str_z(l).paths.ConnDoors(p))-y(l)));


    if y(str_z(l).paths.ConnDoors(p)) -y(l) >...
            str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        a2 =y(str_z(l).paths.ConnDoors(p))...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        b2 =y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        
    elseif y(str_z(l).paths.ConnDoors(p)) -y(l) <...
            -str_z(l).n_bar-str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar)); 
        a2 =-y(l) -str_z(l).n_bar;
        b2 =-y(l) +str_z(l).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    else % |y(str_z(l).paths.ConnDoors(b))-y(l)|<str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(b)).n_bar
     
        a1 =0;
        b1 =sqrt(max(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar),abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar)));
        a2 =min([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        b2 =max([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    end

    
  %%%%% computation of interval with local sensor fault
   fmax=4; fmin=-4;
   f0=0.5*(fmax+fmin);
   fupperbar=0.5*(fmax-fmin);
   flowerbar=-0.5*(fmax-fmin); 
  
  
    if y(str_z(l).paths.ConnDoors(p)) -y(l) + f0 >...
            str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(p)).n_bar + fupperbar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l)+f0 -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar-fupperbar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l)+f0 +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar + fupperbar ));
        a2 =y(str_z(l).paths.ConnDoors(p))...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar + f0;
        b2 =y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar + fupperbar;
        aF =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        bF =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        
    elseif y(str_z(l).paths.ConnDoors(p)) -y(l)+ f0 <...
            -str_z(l).n_bar-str_z(str_z(l).paths.ConnDoors(p)).n_bar+ fupperbar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l)+f0 +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar + fupperbar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l)+f0 -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar - fupperbar)); 
        a2 =-y(l) -str_z(l).n_bar;
        b2 =-y(l) +str_z(l).n_bar;
        aF =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        bF =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    else % |y(str_z(l).paths.ConnDoors(b))-y(l)|<str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(b)).n_bar
     
        a1 =0;
        b1 =sqrt(max(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) + f0 -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar - fupperbar),abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) + f0 +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar + fupperbar )));
        a2 =min([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar )...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar)]);
        b2 =max([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + f0 + fupperbar)]);
        aF =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        bF =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    end

    
    %%%%%
    
    
        Dmu =max(abs(a -c ),abs(b -c ));
        DmuF=max(abs(aF -c),abs(bF-c));    
        sum7 = sum7 + str_z(l).paths.Ad_ij(p) * Dmu;  % h_bar(l);
        
        sum8= sum8 + str_z(l).paths.Ad_ij(p) * DmuF;
        end
        
        
        e1= ((p_air * Cp)/str_z(l).cz) *...
            sqrt( 2 * abs(Cp - Cv)) * sum7;
        e2= ((str_z(l).Umax * awz) / str_z(l).cz) * ...
            (str_z(l).n_bar + n_bar_s) * abs(U(l))...
            +((p_air * Cp)/str_z(l).cz) *...
            sqrt( 2 * abs(Cp - Cv)) * sum7;
        
       
    end
%     if t==0 && flag==0
%         sum88(1,l)=sum7;
%         flag=1;
%     end
    funi1=@(t) exp((A(l,l)- str_z(l).k)* t)*e1;
    funi2=@(t) exp((A(l,l)- str_z(l).k)* t)*e2;
    
    %     sum1=0;
    %     for i=1:tt
    %         sum1=sum1 + str_z(l).lambda^(tt-i)*e1{i};
    %     end
    
    % threshold of zone i
    ebar1(l) = exp((A(l,l)- str_z(l).k)*t)*abs(str_z(l).Tz-str_z(l).Tref)...
            -(1-exp((A(l,l)- str_z(l).k)*t))*(A(l,l)- str_z(l).k)^(-1)...
            *( (abs(str_z(l).k) + ((str_z(l).Umax * awz) / str_z(l).cz) )...
            *str_z(l).n_bar + (1/str_z(l).cz) * sum6 ...
            + ((str_z(l).Umax * awz) / str_z(l).cz)*n_bar_s + str_z(l).r_bar   )...
            + integral(funi1,0,t);
    % + str_z(l).p*sum1;
    % +integral(funi,0,tt);
    %              +lsim(TF1,e1,t);
    %,'RelTol',1e-4,'AbsTol',1e-4);
    
    ebar2(l) = exp((A(l,l)- str_z(l).k)*t)*(abs(str_z(l).Tz-str_z(l).Tref)+0.001)...
            -(1-exp((A(l,l)- str_z(l).k)*t))*(A(l,l)- str_z(l).k)^(-1)...
            *( (abs(str_z(l).k))*str_z(l).n_bar + (1/str_z(l).cz) * sum6 ...
             + str_z(l).r_bar ) + integral(funi2,0,t); 
         
%          ebar2(l)=fval(l);
  
  sigma(l)=(str_z(l).Umax * awz) / str_z(l).cz;
  
   funi3=@(t) exp((A(l,l)- str_z(l).k)* t)*abs(U(l));
%    funi4=@(t) exp((A(l,l)- str_z(l).k)* t)*abs(U(l)-U(l));
%    funi5=@(t) exp((A(l,l)- str_z(l).k)* t)*abs(sum7-sum8);
   funi6=@(t) exp((A(l,l)- str_z(l).k)* t)*sum8;
%    zF(l)= exp((A(l,l)- str_z(l).k)*(t-TD(l)))*abs((str_z(l).k*f0)/(A(l,l)- str_z(l).k))...
%     + ((str_z(l).k*fupperbar)/(A(l,l)- str_z(l).k))*(1-exp((A(l,l)- str_z(l).k)*(t-TD(l))))...
%     + sigma(l)*(abs(f0)+fupperbar)*integral(funi3,TD(l),t)...
%     + sigma(l)*(str_z(l).n_bar + n_bar_s)*integral(funi4,TD(l),t)...
%     + ((p_air * Cp)/str_z(l).cz)*sqrt( 2 * abs(Cp - Cv))*integral(funi5,TD(l),t);

%    zF(l)= exp((A(l,l)- str_z(l).k)*(t-TD(l)))*abs((str_z(l).k*f0)/(A(l,l)- str_z(l).k))...
%     + ( ( str_z(l).k*(fupperbar+ str_z(l).n_bar) + str_z(l).r_bar+ (1/str_z(l).cz) * sum6 )/(A(l,l)- str_z(l).k) )*(1-exp((A(l,l)- str_z(l).k)*(t-TD(l))))...
%     + sigma(l)*(abs(f0) +fupperbar +n_bar_s + str_z(l).n_bar)*integral(funi3,TD(l),t)...
%     + ((p_air * Cp)/str_z(l).cz)*sqrt( 2 * abs(Cp - Cv))*integral(funi6,TD(l),t);



 if DD(l) == 1  
ebarF=0.54;
zF(l)= exp((A(l,l)- str_z(l).k)*(t))*abs(ebarF+(str_z(l).k*f0)/(A(l,l)- str_z(l).k))...
    -( ( str_z(l).k*(fupperbar+str_z(l).n_bar) + str_z(l).r_bar - (1/str_z(l).cz) * sum6 )/(A(l,l)- str_z(l).k) )*(1-exp((A(l,l)- str_z(l).k)*(t)))...
    + sigma(l)*(abs(f0) +fupperbar +n_bar_s + str_z(l).n_bar)*integral(funi3,0,t)...
    +((p_air * Cp)/str_z(l).cz)*sqrt( 2 * abs(Cp - Cv))*integral(funi6,0,t);
eF(l)=(str_z(l).k*abs(f0))/(A(l,l)- str_z(l).k) + zF(l); %+ ebar2(l);

 else
    eF(l)= ebar2(l);    
 end
    
end

else
 eF=zeros(1,Leng);
ebar1=zeros(1,Leng);
ebar2=zeros(1,Leng);
end



if DET==0 
%%
clear i m c
% Adaptive Estimation Scheme FOR ACTUATOR FAULTS
sum1=0;


for i=2*Leng+3:3*Leng+2
    
    m=i-(2*Leng+2);
    if t>= TD(m)
        g(i,i)=((str_z(m).Umax * awz) /(str_z(m).cz)) *...
            (y(Leng+1) - y(m));

        sum1=sum1 + str_z(m).Umax *...
            (y(Leng+1)-y(m))*U(m); %U(m);
        itta(i)=(str_z(m).az/str_z(m).cz)*Ta...
            - ((h*str_z(m).Ad)/str_z(m).cz)*str_z(m).T1;
        sum4=0;
        sum2=0;
        sum3=0;
        for c=1:(length(str_z(m).connected))     %sum for A
            sum4=sum4+str_z(m).az_ij(c)/str_z(m).cz;
            sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)...
                *y(str_z(m).connected(c));
            for d=1:str_z(m).paths.TotalPaths
                sum3 = sum3 +...
                    sign( y(str_z(m).paths.ConnDoors(d))-y(m))*...
                    str_z(m).paths.Ad_ij(d) *...
                    max(y(m),y(str_z(m).paths.ConnDoors(d))) *...
                    sqrt(2 * (Cp - Cv) * ...
                    abs(y(str_z(m).paths.ConnDoors(d)) - y(m)));
                % str_z(m).paths.Ad_ij(d) * x(i) *...
            end
        end
        A(i,i)=(h*str_z(m).Ad-(str_z(m).az))/str_z(m).cz - sum4;
        hh(i) = sum2 + ((p_air * Cp)/str_z(m).cz) * sum3;
        
    else
        g(i,i)=0;
        itta(i)=0;
        A(i,i)=0; 
        hh(i) =0;
 
    end
   Ga(m)=str_z(m).gamma;
   La(m,m)=str_z(m).La;
        
end

% storage tank observer
% A(5*Leng+5, 5*Leng+5)= -aw/Cw;%aw
% g(5*Leng+5, 5*Leng+5)= (Ustmax/Cw)* (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));
% hh(5*Leng+5) = (awz/Cw)* sum1; %aw
% itta(5*Leng+5) = (aw/Cw) * Tpl;%aw

% storage tank observer gain
%L(Leng+1,Leng+1)=L_s;

% Estimator (zones & storage tank)


LEsta=2*Leng+3:3*Leng+2;
Lomega=3*Leng+3:4*Leng+2; %%5*Leng+6:6*(Leng)+5;
LEa=4*Leng+3:5*(Leng)+2; %%6*Leng+7:7*(Leng)+6;
%LEs=7*Leng+3:8*(Leng)+2;
%G(Leng+1)=gamma_s;
Ea(1:Leng)=y(1:Leng)-(x(LEsta))';

dfhata=Ga*diag(x(Lomega))*diag(Ea)*diag(DD(1:Leng));
% - x(LEa)'
dxea= (x(LEsta)'*A(LEsta,LEsta)...
    + (U(1:Leng)+x(LEa)')*g(LEsta,LEsta)...
    + hh(LEsta) + itta(LEsta)...
    + (y(1:Leng)-x(LEsta)')*La...
    + (dfhata)*diag(x(Lomega)))*diag(DD(1:Leng));
%+ (U(1:Leng)+ x(LE)')*g(LEsta,LEsta)...
domegaa=(x(Lomega)'*(A(LEsta,LEsta)-La)...
   +diag(g(LEsta,LEsta))' )*diag(DD(1:Leng));

   funia=@(t) str_z(l).delta_a*exp(-str_z(l).xia * t)*e1a;

    Ebara(l) = str_z(l).delta_a * exp(-str_z(l).xia * (t-TD(l)))...
        *str_z(l).x_bar + abs(x(Lomega(l)))...
        *(str_z(l).dfa_bar)...
        +str_z(l).n_bar +integral(funia,TD(l),t);

    if abs(Ea(l))>Ebara(l)
        Da(l)=1;
        TDa(l)=min(TDa(l),t);

    else
        Da(l)=0;
    end
    
    if t>= TDa(l)
        DDa(l)=1;
    end
    
    
else
  dxea=zeros(1,Leng);
  domegaa=zeros(1,Leng);
  dfhata=zeros(1,Leng); 
  Ea=zeros(1,Leng);
  Ebara=zeros(1,Leng);
end

%%

if DET==0

clear i m c
% Adaptive Estimation Scheme for SENSOR FAULTS
sum1=0;
for i=5*Leng+3:6*Leng+2
    
    m=i-(5*Leng+2);
    if t>= TD(m)
        g(i,i)=(((str_z(m).Umax * awz)) /(str_z(m).cz)) *...
            (y(Leng+1) - y(m)+x(7*Leng+2+m) );%+ 

        sum1=sum1 + str_z(m).Umax *...
            (y(Leng+1) - y(m) + x(7*Leng+2+m) )*U(m); %
        itta(i)=(str_z(m).az/str_z(m).cz)*Ta...
            - ((h*str_z(m).Ad)/str_z(m).cz)*str_z(m).T1;
        sum4=0;
        sum2=0;
        sum3=0;
        for c=1:(length(str_z(m).connected))     %sum for A
            sum4=sum4+str_z(m).az_ij(c)/str_z(m).cz;
            sum2=sum2+(str_z(m).az_ij(c)/str_z(m).cz)...
                *y(str_z(m).connected(c));
            for d=1:str_z(m).paths.TotalPaths
                sum3 = sum3 +...
                    sign( y(str_z(m).paths.ConnDoors(d))-y(m) + x(7*Leng+2+m))*...
                    str_z(m).paths.Ad_ij(d) *...
                    max(y(m)- x(7*Leng+2+m) ,y(str_z(m).paths.ConnDoors(d))) *...
                    sqrt(2 * (Cp - Cv) * ...
                    abs(y(str_z(m).paths.ConnDoors(d)) - y(m) + x(7*Leng+2+m)));
                % str_z(m).paths.Ad_ij(d) * x(i) *...
            end
        end
        A(i,i)=(h*str_z(m).Ad-(str_z(m).az))/str_z(m).cz - sum4;
        hh(i) = sum2 + ((p_air * Cp)/str_z(m).cz) * sum3;
    else
        g(i,i)=0;
        itta(i)=0;
        A(i,i)=0; 
        hh(i) =0; 
    end
    
    Gs(m)=str_z(m).gammas;
    Ls(m,m)=str_z(m).Ls; 
    sigma(m)=(str_z(m).Umax * awz)...
    /(str_z(m).cz);
end

% storage tank observer
% A(5*Leng+5, 5*Leng+5)= -aw/Cw;%aw
% g(5*Leng+5, 5*Leng+5)= (Ustmax/Cw)* (1+(COPmax-1)*(1-((x(Leng+1)-To))/DTmax));
% hh(5*Leng+5) = (awz/Cw)* sum1; %aw
% itta(5*Leng+5) = (aw/Cw) * Tpl;%aw

% storage tank observer gain
%L(Leng+1,Leng+1)=L_s;

% Estimator (zones & storage tank)

LEsts=5*Leng+3:6*Leng+2;
Lomegas=6*Leng+3:7*Leng+2; %%5*Leng+6:6*(Leng)+5;
LEs=7*Leng+3:8*(Leng)+2; %%6*Leng+7:7*(Leng)+6;

%G(Leng+1)=gamma_s;
Es(1:Leng)=y(1:Leng)-x(LEsts)'-x(LEs)';

dfhats=Gs*diag(x(Lomegas)+eye(Leng,1))*diag(Es)*diag(DD(1:Leng)); %

dxes= (x(LEsts)'*A(LEsts,LEsts)...
    + (U(1:Leng))*g(LEsts,LEsts)...
    + hh(LEsts) + itta(LEsts)...
    + (y(1:Leng)-x(LEsts)'-x(LEs)')*Ls...
    + dfhats*diag(x(Lomegas)))*diag(DD(1:Leng));

% domegas=(x(Lomegas)'*(diag(A(Lomegas-(6*Leng+2)))-Ls)...
%     -diag(Ls)'+U(1:Leng)*diag(sigma))*diag(DD(1:Leng));
domegas=(x(Lomegas)'*(diag(A(LEsts))-Ls)...
    + diag(Ls)'+U(1:Leng)*diag(sigma))*diag(DD(1:Leng));

%%
% Isolation threshold
sum5=0;
for l=1:Leng
    sum5=sum5 + str_z(l).Umax * (n_bar_s + str_z(l).n_bar)...
        * abs(U(l));
    sum6=0;
    sum8=0;
    
    for d=1:(length(str_z(l).connected))     %sum for A
        
        sum6 = sum6 + str_z(l).az_ij(d)...
            * (str_z(l).n_bar  + str_z(d).n_bar);
        sum7=0;
        
        for p=1:str_z(l).paths.TotalPaths
             c=sign(y(str_z(l).paths.ConnDoors(p))-y(l))...
                .*max(y(str_z(l).paths.ConnDoors(p)),y(l))...
                .*sqrt(abs(y(str_z(l).paths.ConnDoors(p))-y(l)));

             cc=sign(y(str_z(l).paths.ConnDoors(p))-y(l)+x(LEs(l)))...
                .*max(y(str_z(l).paths.ConnDoors(p)),y(l)-x(LEs(l)))...
                .*sqrt(abs(y(str_z(l).paths.ConnDoors(p))-y(l)+x(LEs(l))));

    if y(str_z(l).paths.ConnDoors(p)) -y(l) >...
            str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
        -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
        -str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        a2 =y(str_z(l).paths.ConnDoors(p))...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        b2 =y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        
    elseif y(str_z(l).paths.ConnDoors(p)) -y(l) <...
            -str_z(l).n_bar-str_z(str_z(l).paths.ConnDoors(p)).n_bar
        
        a1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar));
        b1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar)); 
        a2 =-y(l) -str_z(l).n_bar;
        b2 =-y(l) +str_z(l).n_bar;
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    else % |y(str_z(l).paths.ConnDoors(b))-y(l)|<str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(b)).n_bar
     
        a1 =0;
        b1 =sqrt(max(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar),abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar)));
        a2 =min([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        b2 =max([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)]);
        a =min([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
        b =max([a1 *a2  a1 *b2  b1 *a2  b1 *b2 ]);
    end
     
        sum8 = sum8 + str_z(l).az_ij(d)...
            * (str_z(l).n_bar  + str_z(d).n_bar);
        sum9=0;

    if y(str_z(l).paths.ConnDoors(p)) -y(l) >...
            str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            + str_z(l).fs_bar
        
        aa1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
        -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
        -str_z(l).n_bar - str_z(l).fs_bar ));
        bb1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar +str_z(l).fs_bar ));
        aa2 =y(str_z(l).paths.ConnDoors(p))...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        bb2 =y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar;
        aa =min([aa1*aa2  aa1*bb2  bb1*aa2  bb1*bb2 ]);
        bb =max([aa1*aa2  aa1*bb2  bb1*aa2  bb1*bb2 ]);
        
    elseif y(str_z(l).paths.ConnDoors(p)) -y(l) <...
            -str_z(l).n_bar ...
            -str_z(str_z(l).paths.ConnDoors(p)).n_bar - str_z(l).fs_bar
        
        aa1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar + str_z(l).fs_bar ));
        bb1 =sqrt(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar - str_z(l).fs_bar )); 
        aa2 =-y(l) -str_z(l).n_bar -str_z(l).fs_bar;
        bb2 =-y(l) +str_z(l).n_bar + str_z(l).fs_bar;
        aa =min([aa1 *aa2  aa1 *bb2  bb1 *aa2  bb1 *bb2 ]);
        bb =max([aa1 *aa2  aa1 *bb2  bb1 *aa2  bb1 *bb2 ]);
    else % |y(str_z(l).paths.ConnDoors(b))-y(l)|<str_z(l).n_bar+str_z(str_z(l).paths.ConnDoors(b)).n_bar
     
        aa1 =0;
        bb1 =sqrt(max(abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) -str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            -str_z(l).n_bar - str_z(l).fs_bar),abs(y(str_z(l).paths.ConnDoors(p))...
            -y(l) +str_z(str_z(l).paths.ConnDoors(p)).n_bar...
            +str_z(l).n_bar + str_z(l).fs_bar)));
        aa2 =min([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar )...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar )...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar )]);
        bb2 =max([-min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar)...
            -max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar)...
            min(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar)...
            max(y(str_z(l).paths.ConnDoors(p))...
            +str_z(str_z(l).paths.ConnDoors(p)).n_bar,y(l)...
            +str_z(l).n_bar + str_z(l).fs_bar)]);
        aa =min([aa1 *aa2  aa1 *bb2  bb1 *aa2  bb1 *bb2 ]);
        bb =max([aa1 *aa2  aa1 *bb2  bb1 *aa2  bb1 *bb2 ]);
    end

        Dmuy =max(abs(aa -cc ),abs(bb -cc ));
            
        sum9 = sum9 + str_z(l).paths.Ad_ij(p) * Dmuy;  % h_bar(l);
        
        Dmu =max(abs(a -c ),abs(b -c ));
            
        sum7 = sum7 + str_z(l).paths.Ad_ij(p) * Dmu;  % h_bar(l);
        end
        
        
        e1a= abs( str_z(l).L ) * str_z(l).n_bar + str_z(l).r_bar...
            + ((str_z(l).Umax * awz) / str_z(l).cz) * ...
            (str_z(l).n_bar + n_bar_s) * (abs(U(l)) + str_z(l).fa_bar)...
            + (1/str_z(l).cz) * sum6...
            + ((p_air * Cp)/str_z(l).cz)...
            *sqrt( 2 * abs(Cp - Cv)) * sum7;
        
        e1s= abs( str_z(l).L ) * str_z(l).n_bar + str_z(l).r_bar...
            + ((str_z(l).Umax * awz) / str_z(l).cz) * ...
            (str_z(l).n_bar + n_bar_s) * abs(U(l)) ...
            + (1/str_z(l).cz) * sum8...
            + ((p_air * Cp)/str_z(l).cz) *...
            sqrt( 2 * abs(Cp - Cv)) * sum9;
        
    end
    
 
%(abs(x(LEa(l)))+ str_z(l).fa_bar)  (str_z(l).dfa_bar)
    funis=@(t) str_z(l).delta_s*exp(-str_z(l).xis * t)*e1s;

    Ebars(l) = str_z(l).delta_s * exp(-str_z(l).xis * (t-TD(l)))...
        *str_z(l).x_bar + (abs(x(Lomegas(l)))+1)...
        *Df(l)...
        +str_z(l).n_bar +integral(funis,TD(l),t);
    % (str_z(l).fs_bar + abs(x(LEs(l))))... (str_z(l).dfs_bar)
%(str_z(l).dfs_bar)
      

%     if t>=TLs(l)
%         y(l)=y(l)-x(LEs(l));
%     end
   


  yV(l)=y(l);  

    if DD(l)==1
        %TD(l)=min(TD(l),t);
        
        yV(l)=str_z(l).Tref;
        
        PP(l)=PP(l)+1;
        if PP(l)<=W   
            H(PP(l),l)=abs(Es(l));
        else
            PP(l)=PP(l)-1;
            H(PP(l)-(W-1):PP(l)-1,l)=H(PP(l)-(W-2):PP(l),l);
            H(PP(l),l)=abs(Es(l));
        end
        [~,~,v]=find(H(:,l));
        Ave(l)=sum(v)/length(find(H(:,l)));
        MAXIM(l)=max(v);

        
        if abs(Es(l)) <= Ebars(l)
            TLs(l)=min(TD(l),t);
        end
    end
        
    if abs(Es(l))>Ebars(l)
        DDs(l)=1;
        TDs(l)=min(TDs(l),t);
    else
        DDs(l)=0;
    end
    
    if t>= TDs(l)
        DDs(l)=1;
        if MAXIM(l)<Ebars(l)%Ave(l)<Ebars(l)
           yv(l)=y(l)-x(LEs(l));
           yV(l)=y(l)-x(LEs(l));
           KK(l)=1;
        else
            KK(l)=0;
            yV(l)=str_z(l).Tref;
        end
    end
    

end


else

  %dxes =zeros(1,Leng);
  domegas=zeros(1,Leng); 
  %dfhats =zeros(1,Leng);
  %Es=zeros(1,Leng);
  Ebars=zeros(1,Leng);

end


% Decision

% for k=1:Leng
%     if abs(Ea(k))>Ebara(k)
%         Da(k)=1;
%         TDa(k)=min(TDa(k),t);
%         
%     else
%         Da(k)=0;
%     end
%     
%     if t>= TDa(k)
%         DDa(k)=1;
%     end
%     if abs(Es(k))>Ebars(k)
%         Ds(k)=1;
%         TDs(k)=min(TDs(k),t);
%         
%     else
%         Ds(k)=0;
%     end
%     
%     if t>= TDs(k)
%         DDs(k)=1;
%     end    
% end

%%


%dx=[dx1 dxe dxea domegaa dfhata dxes domegas dfhats]';

dx=[dx1 dxe dxes dfhats domegaa dfhata dxea domegas]';

%mm=mm+1;
assignin('base','Uvar',U);
evalin('base','Uout(end+1,:)=Uvar;');
assignin('base','eFvar',eF);
evalin('base','eFout(end+1,:)=eFvar;');
assignin('base','ebarvar1',ebar1);
evalin('base','ebarout1(end+1,:)=ebarvar1;');
assignin('base','ebarvar2',ebar2);
evalin('base','ebarout2(end+1,:)=ebarvar2;');
assignin('base','yvar',y);
evalin('base','yout(end+1,:)=yvar;');
assignin('base','yVvar',yV);
evalin('base','yVout(end+1,:)=yVvar;');
assignin('base','fvar',f);
evalin('base','fout(end+1,:)=fvar;');

%assignin('base','fhat',x(LEs)');
%evalin('base','fhatout(end+1,:)=fhat;');
assignin('base','fhat',fhat');
evalin('base','fhatout(end+1,:)=fhat;');

assignin('base','z',z');
evalin('base','zout(end+1,:)=z;');
assignin('base','s',s');
evalin('base','sout(end+1,:)=s;');

assignin('base','favar',fa);
evalin('base','faout(end+1,:)=favar;');
assignin('base','Evar',E);
evalin('base','Eout(end+1,:)=Evar;');
assignin('base','Ebarvar',Ebar);
evalin('base','Ebarout(end+1,:)=Ebarvar;');
assignin('base','Dvar',DD);
evalin('base','Dout(end+1,:)=Dvar;');
assignin('base','Kvar',KK);
evalin('base','Kout(end+1,:)=Kvar;');
assignin('base','Dvara',DDa);
evalin('base','Douta(end+1,:)=Dvara;');
assignin('base','Dvars',DDs);
evalin('base','Douts(end+1,:)=Dvars;');
assignin('base','Evara',Ea);
evalin('base','Eouta(end+1,:)=Evara;');
assignin('base','Evarbara',Ebara);
evalin('base','Ebarouta(end+1,:)=Evarbara;');
assignin('base','Evars',Es);
evalin('base','Eouts(end+1,:)=Evars;');
assignin('base','Evarbars',Ebars);
evalin('base','Ebarouts(end+1,:)=Evarbars;');
assignin('base','T',t);
evalin('base','Tout(end+1)=T;');
assignin('base','rho',rho);
evalin('base','rout(end+1,:)=rho;');

assignin('base','Dg',Dg);
evalin('base','Dgout(end+1,:)=Dg;');
assignin('base','gbar',gbar);
evalin('base','gbarout(end+1,:)=gbar;');

assignin('base','gu1',gu);
evalin('base','guout(end+1,:)=gu1;');
assignin('base','gubar1',gubar);
evalin('base','gubarout(end+1,:)=gubar1;');

end
