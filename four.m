%�����ĸ����⣬australian��SVM���⣬��Ҫ�����ļ�myACT.m��extract.m��accu.m
clear;
%�������ݼ�
load('australian.mat');
tic;
%����extract������ͨ�����ݼ���ö��ι滮�����������ȱ�����4��Ϊ�˷���ѡȡ��ʼ��
[G,h,A,b,Ai,bi,x0]=extract(australian_Tr,australian_Tr_Lb,australian_para_C,4);
%���ݵõ��ľ���ȱ��������������ü�������alpha��alpha�Ƕ�ż�����еĲ���
[alp,lam]=myACT(G,h,A,b,Ai,bi,x0,0,1,0);
j=1;maxj=1;
%n��x�ķ�����Ŀ
[~,n]=size(australian_Tr);

%���ݶ�ż���⹫ʽ��w
w=zeros(n,1);
for i=1:length(alp)
    w=w+(alp(i)*australian_Tr_Lb(i)*australian_Tr(i,:))';
end

accuracy=0;
%���ݶ�ż���⹫ʽ��bb����Ϊ���bb��Ҫһ����0<alpha(j)<c
%���ǵ�Ȼ��ѡ����ʹ��ȷ����ߵ�alpha(j)
%ͨ��ѭ�������ŵ�alpha(j),ÿ�ζ���ѵ������׼ȷ�ȣ�ѡ�����ֵ
for count=1:length(alp)
    if(alp(count)>(1.0e-6) && alp(count)<australian_para_C-(1.0e-6))
        j=count;
    end

bb=australian_Tr_Lb(j);
for i=1:length(alp)
    bb=bb-(alp(i)*australian_Tr_Lb(i)* (australian_Tr(i,:)*australian_Tr(j,:)') );
end
%����accu������׼ȷ��
temp=accu(w,australian_Tr,bb,australian_Tr_Lb);
if(temp > accuracy)
    accuracy=temp;maxj=j;end;
end
%����bb����ѵ���Ͳ��Լ���׼ȷ��
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


