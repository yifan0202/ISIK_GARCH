function p_value = uct(EFR,alpha,T)
% Input: EFR, alpha, T
% Output: p-value of the unconditional coverage test
    n1 = T*EFR;
    LR_uc = 2*log((EFR^n1*(1-EFR)^(T-n1))/(alpha^n1*(1-alpha)^(T-n1)));
    p_value = 1-chi2cdf(LR_uc,1);
end