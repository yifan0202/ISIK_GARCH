function [para, ht, LL, se,htplus] = GARCH(R)
% Input: return data
% Output: para(parameters estimated: ------mu omega, alpha1, beta1)
%              ht(filtered variance); LL(log likelihood); 
%              se(standard error); htplus(predicted ht)

    % Initial values
    mu0 = 0;
    omega0 = 8.8081e-7;
    alpha0 = 0.0682;
    beta0 = 0.9292;
%     mu0 = 5.7087e-4;
%     omega0 = 8.5443e-7;
%     alpha0 = 0.0682;
%     beta0 = 0.9294;
    para0 = [mu0,omega0,alpha0,beta0];
    
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    %options = optimset('Display', 'off','MaxFunEval',15000);
    [para,LL] = fminsearch('ll_Function',para0,options,R);
    LL = -1*LL;
    
    VCV = robustvcv('ll_Function_vcv',para,0,R);
    se = sqrt(diag(VCV));
    [~,ht] = ll_Function(para,R);
   htplus = para(2)+para(3)*(R(end)-para(1))^2+para(4)*ht(end);
end