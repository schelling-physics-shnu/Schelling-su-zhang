clear
Q=20;
N=100;
Pi=0;
p=0; 
T=[0:0.1:1,2:2:10,20,30,40,50];
n=20; 

UMPip=zeros(length(T),length(Pi)); 
for r_k=1:length(T)
    UMp=zeros(n,length(Pi));
    for r_n=1:n
        for r_j=1:length(Pi)
            for r_i=1:length(p)
            UMp(r_n,r_j) = f_starlike2completeconn( Pi(r_j),p(r_i),T(r_k),Q,N );
%             r_i,r_n,r_j,r_k
            end
        end
    end
    % UMpn(r_k,:)=UMp;
    UMPip(r_k,:)=sum(UMp)/n 
end

% UMPip=zeros(length(T),length(p)); 
% for r_k=1:length(T)
%     UMp=zeros(n,length(p));
%     for r_n=1:n
%         for r_j=1:length(p)
%             for r_i=1:length(Pi)
%             UMp(r_n,r_j) = f_starlike2completeconn(Pi(r_i),p(r_j),T(r_k),Q,N);
% %             r_i,r_n,r_j,r_k;
%             end
%         end
%     end
% %     UMpn(r_k,:)=UMp;
%     UMPip(r_k,:)=sum(UMp)/n 
% end



% UMPip=zeros(length(Pi),length(p)); 
% for r_i=1:length(Pi)
%     UMp=zeros(n,length(p));
%     for r_n=1:n
%         for r_j=1:length(p)
%             UMp(r_n,r_j) = f_starlike2completeconn( Pi(r_i),p(r_j),Q,N );
%             r_i,r_n,r_j
%         end
%     end
%     UMPip(r_i,:)=sum(UMp)/n 
% end
