%�����ü���������������Ŀ��ͣ���ʽ�Ͳ���ʽ��Լ��������ϵ������x0�ǳ�ʼ��
%cc��Ŀ�꺯���ĳ����chooseѡ1ʱ���������еĵ�ʽԼ����lagrange����
%chooseѡ2ʱ����ռ䷽��
%���x��lam��val��iter�ֱ������Ž⣬lambda������ֵ�͵�������
function [x,lam,val,iter]=myACT(G,h,A,b,Ai,bi,x0,cc,choose,myshow)
%�������ü�������һ��Լ�����ι滮����:
%min f(x)=0.5*x'*G*x+h'*x+cc, cc�ǳ�����
%s.t.  a_i'*x-b_i=0,(i=1,...,l),
%      a_i'*x-b_i>=0,(i=l+1,...,m)
%chooseѡ�������⣨��ʽԼ�����⣩�Ľⷨ��1ѡ��lagrange����2ѡ����ռ䷨
%���x�����Ž⣬lambda����Ӧ����������
%�����Сֵval��f(x)����������iter����Ϣ
tic;
epsilon=1.0e-6;
k=0; 
x=x0;
%����������
kmax=500;
val=0.5*x'*G*x+h'*x+cc;
%ne�ǵ�ʽԼ������Ŀ��ni�ǲ���ʽԼ������Ŀ
ne=length(b); ni=length(bi); 
%lambda��ʼ��
lam=zeros(ne+ni,1);
%act�������ü�
act=ones(ni,1);
%���������ü�
for i=1:ni
    if ( Ai(:,i)'*x-bi(i) > epsilon ), act(i)=0; end
end
%չʾ��ʼֵ
if(myshow==1)
    fprintf('initial values: ');display(act);display(x0);display(val);end
%%%%%%%%%%%%%
%������
while(k<=kmax)
    %���������
    ACT=A;
    %���������ü�act������ʽԼ������Ai������в��뵽�¾���ACT
    for j=1:ni
        if(act(j)>0), ACT=[ACT, Ai(:,j)]; end
    end
    gk=G*x+h;
    %m���¾����Լ����Ŀ��������ʽ�벻��ʽ�ģ�
    [~,m]=size(ACT);
    %���������⣬�����Թ�cholesky�ֽ����ַ����������
    %�˴���b��Ϊ0����
    if(choose==2)
    [dk, lam]=ZSM(G,gk,ACT,zeros(m,1));end
    if(choose==1)
    [dk, lam]=mylagr(G,gk,ACT,zeros(m,1));end 
    if(choose==3)
    [dk, lam]=haha(G,gk,ACT,zeros(m,1));end 
    %չʾlambda
    if(myshow==1)
        fprintf('k= %d, ',k);display(lam);display(dk);
    end
    %%%%%%%%%%%%%%%
    if(norm(dk)<=epsilon)
        %dkΪ0������
        y=0.0;
        %�ҳ���С��lambdak
        if(length(lam)>ne)
            [y,test]=min(lam(ne+1:length(lam)));
        end
        %lambda(k)(i)>=0,act(i)==1ʱֹͣ����
        if(y>=0)
            judge=0;
        else
            judge=1;
            %�����ü�ȥ��q
            for i=1:ni
                if(act(i)&&(ne+sum(act(1:i)))==test)
                    act(i)=0; break;
                end
            end
            %չʾact
            if(myshow==1)
                fprintf('k= %d, ',k);display(act);end
            %%%%%%%%%%%%%%%%%%%%%
        end
        k=k+1;
    else
        judge=1;
        %�󲽳�
        small=1.0;
        for i=1:ni
            %�ɹ�ʽ�󲽳�alpha
            if( (act(i)==0) && (Ai(:,i)'*dk<0) );
                tm=( bi(i)-Ai(:,i)'*x ) / ( Ai(:,i)'*dk );
                if(tm<small)
                    small=tm; ti=i;
                end
            end
        end
        alp=small;
        x=x+alp*dk;
        %չʾx
        val=0.5*x'*G*x+h'*x+cc;
        if(myshow==1)
                fprintf('k= %d, ',k);display(x);display(val);
        end
        %%%%%%%%%%%%%%%
        if(small<1), act(ti)=1; end
        %չʾact
        if(myshow==1)
                fprintf('k= %d, ',k);display(act);end
        %%%%%%%%%%%%%%%
    end
    if(judge==0), break; end
    k=k+1;
end
%min f(x)��ֵ
val=0.5*x'*G*x+h'*x+cc;
%��������
iter=k;
toc;
end

%���������,����α�����pinv()��ֹ�������,lagrange��
function [x,lambda] = mylagr(G,h,A,b)
oldG=pinv(G);
[~,n]=size(A);
if(n>0)
    rb=A'*oldG*h+b;
    lambda=pinv(A'*oldG*A)*rb;
     x=oldG*(A*lambda-h);
else
    x=-oldG*h;
    lambda=zeros(length(b),1);
end
end

%��������������һ�ַ�ʽ����ռ䷨
function [x,lamda]=ZSM(G,h,A,b)
%con��Լ����Ŀ
con=size(b);
%n��x�ķ�����Ŀ
n=size(h);
%QR�ֽ�
[Q,RR]=qr(A);
Q1=[];
Q2=[];
R=[];
for i=1:con
    Q1=[Q1,Q(:,i)];
    R=[R;RR(i,:)];
end
for i=con+1:n
    Q2=[Q2,Q(:,i)];
end
x0=Q1*inv(R')*b;
GGG=Q2'*G*Q2;
hhh=Q2'*(h+G*x0);
%�����Լ�����Ż�������
[d]=mylagr(GGG,hhh,[],[]);
%���x
x=x0+Q2*d;
%���lambda
lamda=(Q1*inv(R'))'*(G*x+h);
end

%������ȥ��
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

