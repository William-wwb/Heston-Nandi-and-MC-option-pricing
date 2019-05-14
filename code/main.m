%Heston-nandi模型参数估计主程序
rf = 0.02859/244;
parameter0 = [1.7686 0.00001 0.0000043859 0.8733 140.5724];
options = optimoptions( @fminunc,...
                         'Algorithm',   'quasi-newton',...
                         'MaxFunEvals', 1E5,...
                         'MaxIter',     1E5,...
                         'TolFun',      1E-6,...
                         'TolX',        1E-6,...
                         'Display',     'iter'                         );
    
[parameter,fval,~,~,~,hess ] = fmincon( @(parameter) log_likelihood(parameter, x, rf),...
                                             parameter0,...
                                             [],[],[],[],...
                                             [0,0,0,0,0],[],...
                                             [],...
                                             options                    );

lambda = parameter(1);
omega = parameter(2);
alpha = parameter(3);
beta = parameter(4);
gam = parameter(5);


