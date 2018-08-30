%����������⣬sonar��SVM���⣬��Ҫ�����ļ�myACT.m��extract.m��accu.m
clear;
%�������ݼ�
load('sonar.mat');
tic;
%����extract������ͨ�����ݼ���ö��ι滮�����������ȱ�����5��Ϊ�˷���ѡȡ��ʼ��
[G,h,A,b,Ai,bi,x0]=extract(sonar_Tr,sonar_Tr_Lb,sonar_para_C,5);
%���ݵõ��ľ���ȱ��������������ü�������alpha��alpha�Ƕ�ż�����еĲ���
[alp,lam]=myACT(G,h,A,b,Ai,bi,x0,0,1,0);
j=1;maxj=1;
%n��x�ķ�����Ŀ
[~,n]=size(sonar_Tr);

%���ݶ�ż���⹫ʽ��w
w=zeros(n,1);
for i=1:length(alp)
    w=w+(alp(i)*sonar_Tr_Lb(i)*sonar_Tr(i,:))';
end
accuracy=0;
%���ݶ�ż���⹫ʽ��bb����Ϊ���bb��Ҫһ����0<alpha(j)<c
%���ǵ�Ȼ��ѡ����ʹ��ȷ����ߵ�alpha(j)
%ͨ��ѭ�������ŵ�alpha(j),ÿ�ζ���ѵ������׼ȷ�ȣ�ѡ�����ֵ
for count=1:length(alp)
    if(alp(count)>(1.0e-10) && alp(count)<sonar_para_C-(1.0e-10))
        j=count;
    end

bb=sonar_Tr_Lb(j);
for i=1:length(alp)
    bb=bb-(alp(i)*sonar_Tr_Lb(i)* (sonar_Tr(i,:)*sonar_Tr(j,:)') );
end
%����accu������׼ȷ��
temp=accu(w,sonar_Tr,bb,sonar_Tr_Lb);
if(temp>accuracy)
    accuracy=temp;maxj=j;end
end
%����bb����ѵ���Ͳ��Լ���׼ȷ��
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


