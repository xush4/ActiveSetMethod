%求解第四个问题，australian的SVM问题，需要依赖文件myACT.m，extract.m和accu.m
clear;
%导入数据集
load('australian.mat');
tic;
%调用extract函数，通过数据集求得二次规划问题所需矩阵等变量，4是为了方便选取初始点
[G,h,A,b,Ai,bi,x0]=extract(australian_Tr,australian_Tr_Lb,australian_para_C,4);
%根据得到的矩阵等变量，调用起作用集方法求alpha，alpha是对偶问题中的参数
[alp,lam]=myACT(G,h,A,b,Ai,bi,x0,0,1,0);
j=1;maxj=1;
%n是x的分量数目
[~,n]=size(australian_Tr);

%根据对偶问题公式求w
w=zeros(n,1);
for i=1:length(alp)
    w=w+(alp(i)*australian_Tr_Lb(i)*australian_Tr(i,:))';
end

accuracy=0;
%根据对偶问题公式求bb，因为求解bb需要一个，0<alpha(j)<c
%我们当然是选择能使正确率最高的alpha(j)
%通过循环找最优的alpha(j),每次都算训练集的准确度，选出最高值
for count=1:length(alp)
    if(alp(count)>(1.0e-6) && alp(count)<australian_para_C-(1.0e-6))
        j=count;
    end

bb=australian_Tr_Lb(j);
for i=1:length(alp)
    bb=bb-(alp(i)*australian_Tr_Lb(i)* (australian_Tr(i,:)*australian_Tr(j,:)') );
end
%调用accu函数算准确度
temp=accu(w,australian_Tr,bb,australian_Tr_Lb);
if(temp > accuracy)
    accuracy=temp;maxj=j;end;
end
%更新bb，算训练和测试集的准确度
bb=australian_Tr_Lb(maxj);
for i=1:length(alp)
    bb=bb-(alp(i)*australian_Tr_Lb(i)* (australian_Tr(i,:)*australian_Tr(maxj,:)') );
end
accuracy_train=accu(w,australian_Tr,bb,australian_Tr_Lb);
accuracy_test=accu(w,australian_Ts,bb,australian_Ts_Lb);
display(w);
display(bb);
display(lam);
display(accuracy_train);
display(accuracy_test);
toc;


