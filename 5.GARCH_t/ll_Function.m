function [LL, ht] = ll_Function(para0, ret)
% Input: para0(initial value for the parameters)
%           ret(returns)
% Output: LL(log likelihood value)
%              ht(filtered conditional variance)
    mu = para0(1);
    omega = para0(2);
    alpha = para0(3);
    beta = para0(4);
    u = para0(5);
    
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
            
    if (alpha+beta>0.9999999)||(min(ht)<=0)||(alpha<=0)||(beta<0.6)||(beta>=1)||(alpha>=1)||((omega+alpha*(ret(end)-mu)^2+beta*ht(end))<0)
        LL = inf;
        LL_G = inf;
    else
            %LL_G = log(normpdf(ret, mu, sqrt(ht)));
        LL_G = log(tpdf((ret- mu)./sqrt(ht),u)./(sqrt(ht)));
        LL = -sum(LL_G);
    end
end