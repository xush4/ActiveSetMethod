%�����������ݹ�����ι滮����ؾ���G��h�ֱ���Ŀ�꺯���Ķ����һ�������
%A,b,Ai,bi�ֱ�����ʽԼ�������Ͳ���ʽԼ��������ϵ���ͳ���
%�����testA��testB��para�ֱ���ѵ������x��y�ͳͷ���ϵ��
function [G,h,A,b,Ai,bi,x0]=extract(testA,testB,para,exer)
%N����������n��x�ķ�������
[N,n]=size(testA);
G=[];
%����G
for i=1:N
    temp=(testA(i,:)'*testB(i))';
    G=[G;temp];
end
G=G*G';
%����h��A��b��Ai��bi
h=-ones(N,1);
b=zeros(N,1);
A=[];
for i=1:N
    A=[A;testB(i)];
end
bi=[zeros(N,1);-para*ones(N,1)];
Ai=[eye(N),-eye(N)];
x0=[];
%x0��Ϊ��ʼֵҪ����Լ���������������������������Լ��������x0
%���ҵ�x0����Լ�����������ϸ����0���ϸ�С�ڳͷ���
for i=1:N
    if(testB(i)>0)
        if(exer==4) %Australian�����x0ֵ
        x0=[x0;0.305];end
        if(exer==5) %Sonar�����x0ֵ
        x0=[x0;0.0089];end
    end
    if(testB(i)<0)
        if(exer==4)  %Australian�����x0ֵ
        x0=[x0;0.247];end
        if(exer==5)  %Sonar�����x0ֵ
        x0=[x0;0.0077];end
    end
end


