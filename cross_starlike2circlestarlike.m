function [Net] = cross_starlike2circlestarlike(Nodes,p)%p is the probability to connecting leaf nodes（连接叶节点概率） 
rand('state',sum(100*clock));
Netl=zeros(Nodes,Nodes);
Netl(1,:)=1;Netl(:,1)=1;
a=zeros(Nodes,Nodes);
for i=2:Nodes
    a(i,i+1:end)=rand(1,Nodes-i);
    a(i+1:end,i)=a(i,i+1:end)';
end
Netl(find(a<p))=1;
Netl=Netl-diag(diag(Netl));
Net=Netl;