function [x,lam]=haha(G,h,A,b)
[n,m]=size(A);
X=mat2cell(G,[m,n-m],[m,n-m]);
Gbb=X(1);Gnb=X(2);Gbn=X(3);Gnn=X(4);
Gbb=cell2mat(Gbb);Gnb=cell2mat(Gnb);Gbn=cell2mat(Gbn);Gnn=cell2mat(Gnn);
XA=mat2cell(A,[m,n-m]);
Ab=XA(1);An=XA(2);
Ab=cell2mat(Ab);An=cell2mat(An);
Xh=mat2cell(h,[m,n-m]);
hb=Xh(1);hn=Xh(2);
hb=cell2mat(hb);hn=cell2mat(hn);

hhn=hn-An*pinv(Ab)*hb+(Gnb-An*pinv(Ab)*Gbb)*pinv(Ab')*b;
GGn=Gnn-Gnb*pinv(Ab')*An'-An*pinv(Ab)*Gbn+An*pinv(Ab)*Gbb*pinv(Ab')*An';
xn=-pinv(GGn)*hhn;
xb=pinv(Ab')*(b-An'*xn);
x=[];
for i=1:m
    x=[x;xb(i)];end
for i=1:n-m
    x=[x;xn(i)];end
lam=pinv(Ab)*(hb+Gbb*xb+Gbn*xn);
end


