%black-Scholes��ʽ�����׼ȷ��ֵ
function [BS_call, BS_put]=optprice_BS(F0,rf,sigma,strike,T)
%�������ر����ֱ𷵻�ŷʽ���Ǻ�ŷʽ�����ļ�ֵ��  
%futures options
%F0:futures price 
%rf:�޷�������
%sigma:�ɼ۲�����
%strike:ִ�м۸�
%T����������

d1=(log(F0./strike)+sigma^2/2*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
%ŷʽ������Ȩ��B-S��ֵ��
BS_call=exp(-rf*T)*(F0.*normcdf(d1)-strike.*normcdf(d2));
%ŷʽ������Ȩ��B-S��ֵ��
BS_put=exp(-rf*T)*(strike.*normcdf(-d2)-F0.*normcdf(-d1));
