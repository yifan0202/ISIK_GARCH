function [LL, ht] = ll_Function(para0, ret)
% Input: para0(initial value for the parameters)
%           ret(returns)
% Output: LL(log likelihood value)
%              ht(filtered conditional variance)
    mu = para0(1);
    omega = para0(2);
    alpha = para0(3);
    beta = para0(4);
    
    T = size(ret,1);
    ht = nan(T,1);
    zt = nan(T,1);
    ht(1) = var(ret);
    r = ret - mu;
    zt(1) = r(1)/sqrt(ht(1));
    for i =2:T
         ht(i) = omega + alpha*r(i-1)^2 + beta*ht(i-1);
         zt(i) = r(i)/sqrt(ht(i));
    end
            
    if (alpha>=1)||(beta>=1)||(beta<0.6)||(min(ht)<=0)||(alpha<=0)||(alpha+beta>0.9999999)||((omega+alpha*(ret(end)-mu)^2+beta*ht(end))<0)
        LL = inf;
        LL_G = inf;
    else
            %LL_G = log(normpdf(ret, mu, sqrt(ht)));
        LL_G = -0.5*(log(2*pi)+log(ht)+zt.^2);
        LL = -sum(LL_G);
    end
end