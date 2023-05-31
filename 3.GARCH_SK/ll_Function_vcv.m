function [sum_llGCE, ll_GCE] = ll_Function_vcv(para, ret)
% Input: para0(initial value for the parameters)
%           ret(returns)
% Output: LL(log likelihood value)
%              ht(filtered conditional variance)

    mu = para(1);
    omega = para(2);
    alpha1 = para(3);
    alpha2 = para(4);
    beta0 = para(5);
    beta1 = para(6);
    delta2 = para(7);
    gamma0 = para(8);
    gamma1 = para(9);
    theta2 = para(10);
    s1 = para(11);
    k1 = para(12);
    
    T = size(ret,1);
    ht = nan(T,1);
    zt = nan(T,1);
    s = nan(T,1);
    k = nan(T,1);
    
    ht(1) = var(ret);
    zt(1) = (ret(1)-mu)/sqrt(ht(1));
    s(1) = s1;
    k(1) = k1;
    
    r = ret-mu;
    for i =2:T
        ht(i) = omega + alpha1*r(i-1)^2 + alpha2*ht(i-1);
        zt(i) = r(i)/sqrt(ht(i));
        s(i) = beta0+beta1*s(i-1)+delta2*zt(i-1);
        k(i) = gamma0 +gamma1*k(i-1)+theta2*abs(zt(i-1));
    end
    if (min(ht)<0)
        sum_llGCE = inf;
        ll_GCE = inf;
    else
        ll_GCE = logGCE(ret,mu,ht,zt,s,k);
        sum_llGCE = sum(ll_GCE);
    end
end