%算准确度的函数，输入w，x和bb，y根据判决函数sign(w*x+bb)与y比较
function [acc]=accu(w,x,bb,y)
%N是样本数目
N=length(y);
acc=0;
for k=1:N
    if(sign(w'*x(k,:)'+bb)==y(k))
        acc=acc+1;end
end
%acc取了%，如acc=80是准确率达80%
acc=100*acc/N;