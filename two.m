%求解第二个问题，需要依赖文件myACT.m
clear;
%目标函数二次项系数矩阵
G=[4 2 2;2 4 0;2 0 2];
%目标函数一次项系数向量
h=[-8;-6;-4];
%约束函数系数矩阵
A=[-1 -1 -2;eye(3)]';
%约束函数常量矩阵
b=[-3;zeros(3,1)];
%初始值
x0=0.5*ones(3,1);
%目标函数的常量
cc=9;
%调用起作用集方法，第二题没有等式约束，用lagrange方法求解子问题
[x,lam,val,iter]=myACT(G,h,[],[],A,b,x0,cc,1,1);
fprintf('最终的结果是：');
display(x);
display(lam);
display(val);
display(iter);
