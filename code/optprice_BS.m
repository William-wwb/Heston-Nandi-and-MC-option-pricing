%black-Scholes公式计算的准确价值
function [BS_call, BS_put]=optprice_BS(F0,rf,sigma,strike,T)
%两个返回变量分别返回欧式看涨和欧式看跌的价值；  
%futures options
%F0:futures price 
%rf:无风险利率
%sigma:股价波动率
%strike:执行价格
%T：到期期限

d1=(log(F0./strike)+sigma^2/2*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
%欧式看涨期权的B-S估值：
BS_call=exp(-rf*T)*(F0.*normcdf(d1)-strike.*normcdf(d2));
%欧式看跌期权的B-S估值：
BS_put=exp(-rf*T)*(strike.*normcdf(-d2)-F0.*normcdf(-d1));
