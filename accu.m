%��׼ȷ�ȵĺ���������w��x��bb��y�����о�����sign(w*x+bb)��y�Ƚ�
function [acc]=accu(w,x,bb,y)
%N��������Ŀ
N=length(y);
acc=0;
for k=1:N
    if(sign(w'*x(k,:)'+bb)==y(k))
        acc=acc+1;end
end
%accȡ��%����acc=80��׼ȷ�ʴ�80%
acc=100*acc/N;