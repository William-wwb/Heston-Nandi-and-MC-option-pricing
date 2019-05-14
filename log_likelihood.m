%极大似然函数
function [res] = log_likelihood(parameter,x,rf)
r = price2ret(x);
h = r;
Z = r;
lambda = parameter(1);
omega = parameter(2);
alpha = parameter(3);
beta = parameter(4);
gam = parameter(5);

%transform : normalize parameters

h(1) = ( omega + alpha )/( 1 - alpha*gam*gam - beta);
Z(1) = ( r(1) - rf - lambda*h(1) ) / sqrt(h(1));
for i=2:length(Z)
    h(i) = omega + alpha .* ( Z(i-1) - gam .* sqrt(h(i-1)) ).^2 + beta .* h(i-1);
    Z(i) = ( r(i) - rf - lambda.*h(i) ) ./ sqrt(h(i));
end
res = -sum(log( normpdf(Z)./sqrt(h) ));
end
