%clear
%p=1;Pi=0;
%T=0;

rand('state',sum(100*clock));
%Q=20;%number of neighborhoods
%N=100;%number of available location in each neighborhoods

Pself=floor(Q*N*0.4*(1-p));%number of selfish particle
Paltr=ceil(Q*N*0.4*p);%number of altruist particle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Net] = cross_starlike2circlestarlike(Q,Pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qoccu=zeros(1,Q);%the number of particles in give neighborhoods
Qutil=zeros(1,Q);%the utility of neighborhoods
U=0;%total utility
lq=zeros(1,Q);%the potential in block q
Lf=0;%effective free energy(or total potential)
Time=5000000;
% Time=100000;%time of simulation
L=zeros(1,Pself+Paltr);% the location(site) of selfish & altruist particles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始随机分布
i=1;
while i<=Pself+Paltr  % initial allocation of all particles
    L(i)=ceil(rand*Q);
    if Qoccu(L(i))<N
        Qoccu(L(i))=Qoccu(L(i))+1;
    else
        continue
    end
    i=i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始均匀分布
% i=1;j=1;
% while i<=Pself+Paltr  % initial allocation of all particles
%     if Qoccu(j)<(Pself+Paltr)/Q
%         L(i)=j;
%         Qoccu(L(i))=Qoccu(L(i))+1;
%     else
%         j=j+1;
%         continue
%     end
%     i=i+1;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iniQoccu=Qoccu;%initial distribution of particles
for i=1:Q %calculate the utility of each neighborhoods 
    if Qoccu(i)/N<0.5
        Qutil(i)=Qoccu(i)/N*2;
    else
        Qutil(i)=2-Qoccu(i)/N*2;
    end
    U=U+Qoccu(i)*Qutil(i);%calculate the total utility
end
iniU=U;%initial U
%%%%%%%%%%%%%%%free energy calculation%%%%%%%%%%%%%%%%%
for i=1:Q %calculate the potential of each block
    if Qoccu(i)/N<=0.5
        snq=1:Qoccu(i);
        lq(i)=sum(snq/N*2);%第i个block的lq
    else
        snq=floor(N/2+1):Qoccu(i);
        lq(i)=sum((1:floor(N/2))/N*2)+sum(2-snq/N*2);
    end
end
Lf=sum(lq);
iniLf=Lf;%initial Lf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0;
a=0;
while t~=Time
%     rp=ceil(rand*(Pself+Paltr));%choose a particle randomly
%     rQ=ceil(rand*Q);%choose a neighborhood randomly
    roL=find(L==ceil(rand*Q));
    if length(roL)==0,continue,end
    rp=roL(ceil(rand*length(roL)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rlink=find(Net(L(rp),:)>0);%rlink is the neighbor of node L(rp)
    rQ=rlink(ceil(rand*length(rlink)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Qoccu(rQ)==N | L(rp)==rQ,continue,end %the chosen neighborhood should not be full
    rQutil=Qutil;%rQutil is a temperary Qutil
    rQoccu=Qoccu;%rQoccu is a temperary Qoccu
    rQoccu(rQ)=rQoccu(rQ)+1; rQoccu(L(rp))=rQoccu(L(rp))-1;
    if rQoccu(rQ)/N<0.5
        rQutil(rQ)=rQoccu(rQ)/N*2;
    else
        rQutil(rQ)=2-rQoccu(rQ)/N*2;
    end
    if rQoccu(L(rp))/N<0.5
        rQutil(L(rp))=rQoccu(L(rp))/N*2;
    else
        rQutil(L(rp))=2-rQoccu(L(rp))/N*2;
    end
    rU=sum(rQutil.*rQoccu);%temperary U
    
    F=rand;
    
    if rp<=Pself %to judge the chosen particle is selfish or the other 
%         
        PR=1/(1+exp(-(rQutil(rQ)-Qutil(L(rp)))/T));
    if F<PR
%         if Qutil(L(rp))<rQutil(rQ)
           
            Qoccu=rQoccu;
            Qutil=rQutil;
            U=rU;
            L(rp)=rQ;
        end
    else
        PR=1/(1+exp(-(rU-U)/T));
        if F<PR
%         if U<rU
     
            Qoccu=rQoccu;
            Qutil=rQutil;
            U=rU;
            L(rp)=rQ;
        end
    end
%     U
    t=t+1;
    if rem(t,1000)==0
        t;
%         Qoccu(2)
        [mQoccu,lQoccu]=max(Qoccu);
        oQoccu=find(Qoccu>0);
        [mQutil,lQutil]=min(Qutil(oQoccu));
        Qoccu(oQoccu(lQutil));
%         lQutil
        Qoccu(lQutil);
        sQoccu(t/1000,:)=Qoccu;
        sQutil(t/1000,:)=Qutil;
    end
    if t>1000000 && rem(t,10000)==0
        a=a+1;
        for i=1:Q %calculate the potential of each block
            if Qoccu(i)/N<=0.5
                snq=1:Qoccu(i);
                lq(i)=sum(snq/N*2);
            else
                snq=floor(N/2+1):Qoccu(i);
                lq(i)=sum((1:floor(N/2))/N*2)+sum(2-snq/N*2);
            end
        end
        Lf=sum(lq);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aum(a,1)=U/(Pself+Paltr);% average utility
        Lfm(a,1)=Lf/(Pself+Paltr);% average effective free energy
    end
    if a>6 && abs(aum(a,1)-aum(a-1,1))<0.005 && abs(aum(a-1,1)-aum(a-2,1))<0.005 ...
            && abs(aum(a-2,1)-aum(a-3,1))<0.005 && abs(aum(a-3,1)-aum(a-4,1))<0.005...
            && abs(aum(a-4,1)-aum(a-5,1))<0.005
        break
    end
end


%%%%%%%%%%%%%%free energy calculation%%%%%%%%%%%%%%%%%%%
for i=1:Q %calculate the potential of each block 
            if Qoccu(i)/N<=0.5
                snq=1:Qoccu(i);
                lq(i)=sum(snq/N*2);
            else
                snq=floor(N/2+1):Qoccu(i);
                lq(i)=sum((1:floor(N/2))/N*2)+sum(2-snq/N*2);
            end
end
t
Lf=sum(lq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
um=U/(Pself+Paltr);% average utility
Lfm=Lf/(Pself+Paltr);% average effective free energy
 
    
