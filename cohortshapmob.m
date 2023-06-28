function [Shap]=cohortshapmob(x,y,d)
% Shapley values via Möbius transform. 
% Original file written by Elmar Plischke ("Computing Shapley Effects for Sensitivity Analysis" with E. Borgonovo and G. Rabitti, SIAM/ASA JUQ, 2021) 
% with modified value function estimation using cohorts.  
%%

k=size(x,2);
n=size(x,1);
l=2^k-1;H=zeros(1,l);
sz=zeros(1,l);
for i=1:l 
 g=bitand(i,2.^(0:k-1))~=0;
 m1=sum(g); 
 sz(i)=m1;
 ybar=zeros(1,n);
 for t=1:n
     C = 0;
     index = abs(x(:,g)-x(t,g))<=d(g);
     C=(sum(index,2)==m1);
     ybar(t)=mean(y(find(C==1)));
 end 
 H(i)=mean((ybar-mean(y)).^2); 
end
mob=zeros(size(H));
sel=false(1,l);
for i=1:l
    sel(1:i)=xor(sel(1:i),[true,sel(1:i-1)]);
    ii=find(sel(1:i));
    mob(:,i)=(H(:,ii)*(-1).^(sz(i)+sz(ii)'))/sz(i);
end
Shapunorm=ones(size(H,1),k);
for i=1:k
    Shapunorm(:,i)=sum(mob(:,logical(bitand(1:l,2^(i-1)))),2);
end
%% variance
Shap=Shapunorm./H(:,end);