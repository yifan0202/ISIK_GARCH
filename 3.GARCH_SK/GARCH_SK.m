function [para, ht, zt,s,k,LL, ll_r,se,htplus,splus,kplus] = GARCH_SK(ret)
% Input: return data
% Output: para(parameters estimated: ------mu omega, alpha1, beta1,)
%              ht(filtered variance); LL(log likelihood); 
%              se(standard error); htplus(predicted ht)
    % initial values

    mu = 5.7087e-4;
    omega = 8.5443e-7;
    alpha1 = 0.0682;
    alpha2 = 0.9294;
    
    beta0 = -0.0160;
    beta1 = -0.2000;
    delta2 = 0.05;
    
    gamma0 = 0.004;
    gamma1 = 0.95;
    theta2 = 0.04;
    
    s1 = 0.1761;
    k1 = 5.5973;
%     mu = 0;
%     omega = 8.8081e-7;
%     alpha1 = 0.0682;
%     alpha2 = 0.9292;
%     
%     beta0 = -0.0160;
%     beta1 = -0.378;
%     delta2 = 0.162;
%     
%     gamma0 = 0.004;
%     gamma1 = 0.822;
%     theta2 = 0.086;
%     
%     s1 = 0.1761;
%     k1 = 5.5973;
    
    para0 = [mu, omega, alpha1, alpha2, ...
        beta0, beta1, delta2,...
        gamma0, gamma1, theta2,...
        s1, k1];
    
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    [para,LL] = fminsearch('ll_Function',para0,options,ret);
    LL = -1*LL;
    
    VCV = robustvcv('ll_Function_vcv',para,0,ret);
    se = sqrt(diag(VCV));
    [~,ll_r,ht,zt,s,k] = ll_Function(para,ret);
   htplus = para(2)+para(3)*(ret(end)-para(1))^2+para(4)*ht(end);
   splus = para(5)+para(6)*s(end)+para(7)*zt(end);
   kplus = para(8)+para(9)*k(end)+para(10)*abs(zt(end));
   
end