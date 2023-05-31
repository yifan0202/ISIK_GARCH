function ll = logGCE(ret,mu,ht,zt,s,k)
% Input: ret, mu, ht, zt, s, k
% Output: loglikelihood function of rt with transformed GCE distribution
    Phi = 1+ s/6.*(zt.^3-3*zt)+(k-3)./24.*(zt.^4-6*zt.^2+3);
    Gamma = 1+s./6+(k-3).^2./24;
    ll = log(normpdf(ret,mu,sqrt(ht)))+log(Phi.^2)-log(Gamma);
end