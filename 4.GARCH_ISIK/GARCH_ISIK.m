function [para, ht,zt,s,k,eta,vartheta,LL,LL_r,LL_IS,LL_IK, se, htplus,splus,kplus] = GARCH_ISIK(ret, ISIK)
% Input:ret(return), ISIK(implied skewness and implied kurtosis)
% Output: para(parameters)
%              ht(conditional variance)
%              zt, s, k, eta, vartheta
%              LL = LL_GCE+LL_IS+LL_IK
%              se: robust standard error
%              htplut, splus,kplus (predicted ht, s and k)
    IS = ISIK(:,1);
    IK = ISIK(:,2);

%     mu = 0;
%     omega = 8.8081e-7;
%     alpha1 = 0.0682;
%     alpha2 = 0.9292;
%     
%     beta0 = -0.0130;
%     beta1 = -0.0100;
%     beta2 = 0.2;
%     
%     delta0 = 0;
%     delta1 = 0;
%     delta2 = 0.065;
%     
%     gamma0 = 0.055;
%     gamma1 = 0.964;
%     gamma2 = 0;
%     
%     theta0 = 0;
%     theta1 = 0;
%     theta2 = -0.03;
%     
%     std_eta = 1;
%     std_vartheta = 0.3410;
%     
%     s1 = 0.1761;
%     k1 = 5.5973;
%     mu = 0;
%     omega = 8.8081e-7;
%     alpha1 = 0.0682;
%     alpha2 = 0.9292;
    mu = 5.7087e-4;
    omega = 8.5443e-7;
    alpha1 = 0.0682;
    alpha2 = 0.9294;
    
    beta0 = -0.0130;
    beta1 = 0.046;
    beta2 = 0.086;
    
    delta0 = 0;
    delta1 = -0.771;
    delta2 = 0.562;
    
    gamma0 = 0.055;
    gamma1 = 0.974;
    gamma2 = 0.009;
    
    theta0 = 0;
    theta1 = 2.336;
    theta2 = 0.060;
    
    std_eta = 1;
    std_vartheta = 0.3410;
    
    s1 = 0.1761;
    k1 = 5.5973;
    
    para0 = [mu,omega,alpha1,alpha2,...
        beta0,beta1,beta2,...
        delta0, delta1, delta2,...
        gamma0,gamma1,gamma2,...
        theta0, theta1,theta2,...
        std_eta,std_vartheta,...
        s1,k1];
    
    data = [ret,IS,IK];
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    [para,LL] = fminsearch('ll_Function',para0,options,data);
    LL = -1*LL;
    
    VCV = robustvcv('ll_Function_vcv',para,0,data);
    se = sqrt(diag(VCV));
    [~,LL_r,LL_IS,LL_IK,ht,zt,s,k,eta,vartheta] = ll_Function(para,data);
    htplus = para(2)+para(3)*(ret(end)-para(1))^2+para(4)*ht(end);
    splus = para(5)+para(6)*s(end)+para(7)*IS(end);
    kplus = para(11)+para(12)*k(end)+para(13)*IK(end);
end