function OptionPrice=HestonNandi(S_0,X,Sig_,T,r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the price of Call option based on the GARCH 
% option pricing formula of Heston and Nandi(2000). The input to the
% function are: current price of the underlying asset, strike price,
% unconditional variance of the underlying asset, time to maturity in days,
% and daily risk free interest rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ali Boloorforoosh
% email:  a_bol@jmsb.concordia.ca
% Date:   Nov. 1,08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%% sample inputs %%%%%
    % S_0=100;                    stock price at time t
    % X=100;                      strike prices
    % Sig_=.04/252;               unconditional variances per day
    % T=30;                       option maturity
    % r=.05/365;                  daily risk free rate


OptionPrice=.5*S_0+(exp(-r*T)/pi)*quad(@Integrand1,eps,100)-X*exp(-r*T)*(.5+(1/pi)*quad(@Integrand2,eps,100));

    % function Integrand1 and Integrand2 return the values inside the 
    % first and the second integrals
    
    function f1=Integrand1(phi)
        f1=real((X.^(-i*phi).*charac_fun(i*phi+1))./(i*phi));
    end

    function f2=Integrand2(phi)
        f2=real((X.^(-i*phi).*charac_fun(i*phi))./(i*phi));
    end

    % function that returns the value for the characteristic function
    function f=charac_fun(phi)
        
        phi=phi';    % the input has to be a row vector
        
        % GARCH parameters
        lam_=-.5;                   % risk neutral version of lambda
    %soybean meal
%         lam=0.018226894135569;
%         a=5.299413743790637e-09;
%         b=0.500234453425116;
%         g=1.437432176713487e+02;
%         w=1.183670812003007e-04;
    %corn
%         lam=0.001347743604854;
%         a=2.527694665888083e-05;
%         b=0.722794516716392;
%         g=9.052885514341369e-05;
%         w=1.297183134837510e-05;

    % white sugar
          lam=1.59587173568312;
          a=8.355814757916855e-09;
          b=0.499963298593962;
          g=1.405578310141568e+02;
%           w=3.873841932272978e-05;
    %cotton No.1
%         lam=1.566986967219778;
%         a=9.906805531178909e-07;
%         b=0.970787024934775;
%         g=1.403933798984022e+02;
%         w=3.287949638107078e-08;
        
          
%         lam=1.7686;
%         a=.0000043859;
%         b=.8733;
%         g=140.5724;                      % gamma coefficient
        g_=g+lam+.5;                % risk neutral version of gamma
        w=Sig_*(1-b-a*g^2)-a;       % GARCH intercept
        
        % recursion for calculating A(t,T,Phi)=A_ and B(t,T,Phi)=B_
        A(:,T-1)=phi.*r;
        B(:,T-1)=lam_.*phi+.5*phi.^2;

        for i=2:T-1
            A(:,T-i)=A(:,T-i+1)+phi.*r+B(:,T-i+1).*w-.5*log(1-2*a.*B(:,T-i+1));
            B(:,T-i)=phi.*(lam_+g_)-.5*g_^2+b.*B(:,T-i+1)+.5.*(phi-g_).^2./(1-2.*a.*B(:,T-i+1));
        end

        A_=A(:,1)+phi.*r+B(:,1).*w-.5*log(1-2.*a.*B(:,1));                    % A(t;T,phi)
        B_=phi.*(lam_+g_)-.5*g_^2+b.*B(:,1)+.5*(phi-g_).^2./(1-2.*a.*B(:,1)); % B(t;T,phi)

        f=S_0.^phi.*exp(A_+B_.*Sig_);
        f=f'; % the output is a row vector

    end

end





