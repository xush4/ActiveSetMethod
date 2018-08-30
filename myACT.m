%起作用集方法函数，输入目标和（等式和不等式）约束函数的系数矩阵，x0是初始点
%cc是目标函数的常数项，choose选1时，子问题中的等式约束用lagrange方法
%choose选2时用零空间方法
%输出x，lam，val和iter分别是最优解，lambda，最优值和迭代次数
function [x,lam,val,iter]=myACT(G,h,A,b,Ai,bi,x0,cc,choose,myshow)
%用起作用集方法解一般约束二次规划问题:
%min f(x)=0.5*x'*G*x+h'*x+cc, cc是常数项
%s.t.  a_i'*x-b_i=0,(i=1,...,l),
%      a_i'*x-b_i>=0,(i=l+1,...,m)
%choose选择子问题（等式约束问题）的解法，1选择lagrange法，2选择零空间法
%输出x是最优解，lambda是相应乘子向量；
%输出极小值val是f(x)，迭代次数iter等信息
tic;
epsilon=1.0e-6;
k=0; 
x=x0;
%最大迭代次数
kmax=500;
val=0.5*x'*G*x+h'*x+cc;
%ne是等式约束的数目，ni是不等式约束的数目
ne=length(b); ni=length(bi); 
%lambda初始化
lam=zeros(ne+ni,1);
%act是起作用集
act=ones(ni,1);
%更新起作用集
for i=1:ni
    if ( Ai(:,i)'*x-bi(i) > epsilon ), act(i)=0; end
end
%展示初始值
if(myshow==1)
    fprintf('initial values: ');display(act);display(x0);display(val);end
%%%%%%%%%%%%%
%迭代层
while(k<=kmax)
    %求解子问题
    ACT=A;
    %根据起作用集act将不等式约束矩阵Ai的相关列插入到新矩阵ACT
    for j=1:ni
        if(act(j)>0), ACT=[ACT, Ai(:,j)]; end
    end
    gk=G*x+h;
    %m是新矩阵的约束数目（包括等式与不等式的）
    [~,m]=size(ACT);
    %子问题的求解，但尝试过cholesky分解会出现非正定的情况
    %此处将b变为0向量
    if(choose==2)
    [dk, lam]=ZSM(G,gk,ACT,zeros(m,1));end
    if(choose==1)
    [dk, lam]=mylagr(G,gk,ACT,zeros(m,1));end 
    if(choose==3)
    [dk, lam]=haha(G,gk,ACT,zeros(m,1));end 
    %展示lambda
    if(myshow==1)
        fprintf('k= %d, ',k);display(lam);display(dk);
    end
    %%%%%%%%%%%%%%%
    if(norm(dk)<=epsilon)
        %dk为0的情形
        y=0.0;
        %找出最小的lambdak
        if(length(lam)>ne)
            [y,test]=min(lam(ne+1:length(lam)));
        end
        %lambda(k)(i)>=0,act(i)==1时停止迭代
        if(y>=0)
            judge=0;
        else
            judge=1;
            %起作用集去掉q
            for i=1:ni
                if(act(i)&&(ne+sum(act(1:i)))==test)
                    act(i)=0; break;
                end
            end
            %展示act
            if(myshow==1)
                fprintf('k= %d, ',k);display(act);end
            %%%%%%%%%%%%%%%%%%%%%
        end
        k=k+1;
    else
        judge=1;
        %求步长
        small=1.0;
        for i=1:ni
            %由公式求步长alpha
            if( (act(i)==0) && (Ai(:,i)'*dk<0) );
                tm=( bi(i)-Ai(:,i)'*x ) / ( Ai(:,i)'*dk );
                if(tm<small)
                    small=tm; ti=i;
                end
            end
        end
        alp=small;
        x=x+alp*dk;
        %展示x
        val=0.5*x'*G*x+h'*x+cc;
        if(myshow==1)
                fprintf('k= %d, ',k);display(x);display(val);
        end
        %%%%%%%%%%%%%%%
        if(small<1), act(ti)=1; end
        %展示act
        if(myshow==1)
                fprintf('k= %d, ',k);display(act);end
        %%%%%%%%%%%%%%%
    end
    if(judge==0), break; end
    k=k+1;
end
%min f(x)的值
val=0.5*x'*G*x+h'*x+cc;
%迭代次数
iter=k;
toc;
end

%求解子问题,采用伪逆矩阵pinv()防止奇异矩阵,lagrange法
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

%求解子问题的另外一种方式，零空间法
function [x,lamda]=ZSM(G,h,A,b)
%con是约束数目
con=size(b);
%n是x的分量数目
n=size(h);
%QR分解
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
%求解无约束最优化子问题
[d]=mylagr(GGG,hhh,[],[]);
%算得x
x=x0+Q2*d;
%算得lambda
lamda=(Q1*inv(R'))'*(G*x+h);
end

%变量消去法
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

