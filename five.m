%求解第五个问题，sonar的SVM问题，需要依赖文件myACT.m，extract.m和accu.m
clear;
%导入数据集
load('sonar.mat');
tic;
%调用extract函数，通过数据集求得二次规划问题所需矩阵等变量，5是为了方便选取初始点
[G,h,A,b,Ai,bi,x0]=extract(sonar_Tr,sonar_Tr_Lb,sonar_para_C,5);
%根据得到的矩阵等变量，调用起作用集方法求alpha，alpha是对偶问题中的参数
[alp,lam]=myACT(G,h,A,b,Ai,bi,x0,0,1,0);
j=1;maxj=1;
%n是x的分量数目
[~,n]=size(sonar_Tr);

%根据对偶问题公式求w
w=zeros(n,1);
for i=1:length(alp)
    w=w+(alp(i)*sonar_Tr_Lb(i)*sonar_Tr(i,:))';
end
accuracy=0;
%根据对偶问题公式求bb，因为求解bb需要一个，0<alpha(j)<c
%我们当然是选择能使正确率最高的alpha(j)
%通过循环找最优的alpha(j),每次都算训练集的准确度，选出最高值
for count=1:length(alp)
    if(alp(count)>(1.0e-10) && alp(count)<sonar_para_C-(1.0e-10))
        j=count;
    end

bb=sonar_Tr_Lb(j);
for i=1:length(alp)
    bb=bb-(alp(i)*sonar_Tr_Lb(i)* (sonar_Tr(i,:)*sonar_Tr(j,:)') );
end
%调用accu函数算准确度
temp=accu(w,sonar_Tr,bb,sonar_Tr_Lb);
if(temp>accuracy)
    accuracy=temp;maxj=j;end
end
%更新bb，算训练和测试集的准确度
bb=sonar_Tr_Lb(maxj);
for i=1:length(alp)
    bb=bb-(alp(i)*sonar_Tr_Lb(i)* (sonar_Tr(i,:)*sonar_Tr(maxj,:)') );
end
accuracy_train=accu(w,sonar_Tr,bb,sonar_Tr_Lb);
accuracy_test=accu(w,sonar_Ts,bb,sonar_Ts_Lb);
display(w);
display(bb);
display(lam);
display(accuracy_train);
display(accuracy_test);
toc;


