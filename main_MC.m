%鸡蛋期货monte-carlo定价主程序
s0 = x(length(x));
rf = 0.04/252;
strike = s0;
T = 100;
M = 100000;
[HN_call,HN_put]=optprice_HN(s0,strike,rf,lambda,alpha,beta,gam,omega,T,M)

