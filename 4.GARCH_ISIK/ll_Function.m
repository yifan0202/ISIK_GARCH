function [sum_ll,LL_r,LL_IS,LL_IK,ht,zt,s,k,eta,vartheta] = ll_Function(para, data)
% Input: para(parameters estimated), ret(ret)
% Output: sum_ll(sum of the log likelihood ratio)
%              LL_IS (log likelihood ratio for the implied skewness)
%              LL_IK (log likelihood ratio for the implied skewness)
%              ht, zt, s, k, eta, vartheta
    ret = data(:,1);
    IS = data(:,2);
    IK = data(:,3);

    mu = para(1);
    omega = para(2);
    alpha1 = para(3);
    alpha2 = para(4);
    beta0 = para(5);
    beta1 = para(6);
    beta2 = para(7);
    delta0 = para(8);
    delta1 = para(9);
    delta2 = para(10);
    gamma0 = para(11);
    gamma1 = para(12);
    gamma2 = para(13);
    theta0 = para(14);
    theta1 = para(15);
    theta2 = para(16);
    std_eta = para(17);
    std_vartheta = para(18);
    s1 = para(19);
    k1 = para(20);
    
    T = size(ret,1);
    ht = nan(T,1);
    zt = nan(T,1);
    s = nan(T,1);
    k = nan(T,1);
    
    ht(1) = var(ret);
    r = ret-mu;%demean
    zt(1) = r(1)/sqrt(ht(1));
    s(1) = s1;
    k(1) = k1;
    
    for i = 2:T
        ht(i) = omega + alpha1*r(i-1)^2 + alpha2*ht(i-1);
        zt(i) = r(i)/sqrt(ht(i));
        s(i) = beta0+beta1*s(i-1)+beta2*IS(i-1);
        k(i) = gamma0 +gamma1*k(i-1)+gamma2*IK(i-1); 
    end
    eta = IS-delta0-delta1*s-delta2*zt;
    vartheta = IK - theta0-theta1*k - theta2*abs(zt);
    [~, u2, ~, u4] = Transf(s,k);
   
    if (min(ht)<=0)||(min(k)<1)||(alpha1+alpha2>1)||(std_eta<0)||(std_vartheta<0)||(sum(u2<0)>0)||(sum(u4<1)>0)
        sum_ll = inf;
    else
        LL_GCE = logGCE(ret,mu,ht,zt,s,k);
        LL_IS = sum(-0.5*(log(2*pi)+log(std_eta.^2)+log(eta.^2./std_eta.^2)));
        LL_IK = sum(-0.5*(log(2*pi)+log(std_vartheta.^2)+log(vartheta.^2./std_vartheta.^2)));
        sum_ll = -(sum(LL_GCE)+LL_IS+LL_IK);
        if sum_ll == -inf
            sum_ll = inf;
        end
    end
    LL_r = sum(-0.5*(log(2*pi)+log(ht)+zt.^2));
end