function [para, ht, LL, se,htplus] = GARCH_t(R)
% Input: return data
% Output: para(parameters estimated: ------mu omega, alpha1, beta1)
%              ht(filtered variance); LL(log likelihood); 
%              se(standard error); htplus(predicted ht)

    % Initial values
%     mu0 = 0;
%     omega0 = 8.8081e-7;
%     alpha0 = 0.0682;
%     beta0 = 0.9292;
    mu0 = 0;
    omega0 = 1.7532e-6;
    alpha0 = 0.0985;
    beta0 = 0.9013;
    u0 = 4.1437;
    para0 = [mu0,omega0,alpha0,beta0,u0];
    
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    %options = optimset('Display', 'off','MaxFunEval',15000);
    [para,LL] = fminsearch('ll_Function',para0,options,R);
    LL = -1*LL;
    
    VCV = robustvcv('ll_Function_vcv',para,0,R);
    se = sqrt(diag(VCV));
    [~,ht] = ll_Function(para,R);
   htplus = para(2)+para(3)*(R(end)-para(1))^2+para(4)*ht(end);
end