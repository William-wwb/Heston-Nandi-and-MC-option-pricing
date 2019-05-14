function[Eur_call,Eur_put]=optprice_HN(s0,strike,rf,lambda,alpha,beta,gam,omega,T,M)
%s0Ϊ��ʼ�۸�
%strikeΪִ�м�
%rΪ�޷�������
%lambda alpha beta gamma omega ��ΪHN���Ƴ����Ĳ���
%T��������
%N����
%Mģ�����
sigma0 = (omega+alpha)/(1-beta-alpha*gam.^2);
gam = gam + lambda + 0.5;
lambda = -1/2;
s = repmat(s0,M,1);
h = repmat(sigma0,M,1);
e = randn(M,T+1);
for i=1:T
    h =[h (omega+beta.*h(:,i)+alpha.*(e(:,i)-gam.*sqrt(h(:,i))).^2)];
    s =[s s(:,i).*exp(rf-1/2*h(:,i+1)+sqrt(h(:,i+1)).*e(:,i+1))];
end
%������Ȩ
cpayoff = exp(-rf*T).*max(s(:,T+1)-strike,0);
Eur_call = normfit(cpayoff,0.01);
%������Ȩ
ppayoff = exp(-rf*T)*max(strike-s(:,T+1),0);
Eur_put = normfit(ppayoff,0.01);
end









    