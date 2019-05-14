%price cotton No.1 options on 26th April 2019 
load('fourFutures.mat');
x=cotton;
%first, run main.m to get parameters
s0=15335;
rf=0.02859;
T=25/244;  %   years
strike=(14000:200:17000)';
N=size(strike,1);
sigma=std(price2ret(x))*sqrt(244);
%%  BS Models
[~, BS_put]=optprice_BS(s0,rf,sigma,strike,T);
%%   solution of HN
T2=27;  %day
r2=rf/365;                  %daily risk free rate
sigma2= (omega+alpha)/(1-beta-alpha*gam.^2);%unconditional variances per day
Call_HN=zeros(N,1);
Put_HN=zeros(N,1);
for i=1:N
    Call_HN(i)=HestonNandi(s0,strike(i),sigma2,T2,r2);
    Put_HN(i)=Call_HN(i)+(strike(i)-s0)*exp(-r2*T2);
end
%%    MC of HN European put options  
T3=27;%day
r3=rf/365;                  %daily risk free rate
M = 100000;
N=size(strike,1);
MC_put=zeros(N,1);
for i=1:N
[~,MC_put(i)]=optprice_HN(s0,strike(i),r3,lambda,alpha,beta,gam,omega,T3,M);
end