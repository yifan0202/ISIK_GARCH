function VaR = getVaR_sk(mu,s,k,alpha,htplus)
% Input: mu, s, k; alpha(level of confidence),htplus
% Output: VaR
    phi_s = norminv(alpha)+(norminv(alpha)^2-1)*s/6+(norminv(alpha)^3-3*norminv(alpha))*(k-3)/24;
    VaR = -(mu+sqrt(htplus)*phi_s);
end