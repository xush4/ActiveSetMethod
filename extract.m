%根据所得数据构造二次规划的相关矩阵，G，h分别是目标函数的二次项，一次项矩阵
%A,b,Ai,bi分别代表等式约束函数和不等式约束函数的系数和常量
%输入的testA，testB和para分别是训练集的x，y和惩罚项系数
function [G,h,A,b,Ai,bi,x0]=extract(testA,testB,para,exer)
%N是样本数，n是x的分量个数
[N,n]=size(testA);
G=[];
%构造G
for i=1:N
    temp=(testA(i,:)'*testB(i))';
    G=[G;temp];
end
G=G*G';
%构造h，A，b，Ai和bi
h=-ones(N,1);
b=zeros(N,1);
A=[];
for i=1:N
    A=[A;testB(i)];
end
bi=[zeros(N,1);-para*ones(N,1)];
Ai=[eye(N),-eye(N)];
x0=[];
%x0作为初始值要满足约束条件，这里根据样本找了满足约束条件的x0
%我找的x0满足约束条件，且严格大于0和严格小于惩罚项
for i=1:N
    if(testB(i)>0)
        if(exer==4) %Australian所需的x0值
        x0=[x0;0.305];end
        if(exer==5) %Sonar所需的x0值
        x0=[x0;0.0089];end
    end
    if(testB(i)<0)
        if(exer==4)  %Australian所需的x0值
        x0=[x0;0.247];end
        if(exer==5)  %Sonar所需的x0值
        x0=[x0;0.0077];end
    end
end


