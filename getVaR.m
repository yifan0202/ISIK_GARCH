function VaR = getVaR(mu, u, alpha,htplus)
% Input: mu, u(u1, u2, u3, u4); alpha(level of confidence),htplus
% Output: VaR
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
    u4 = u(4);
    
    u3_s = u3/(u2^(3/2));
    u4_s = u4/(u2^2);
    
    phi_s = norminv(alpha)+(norminv(alpha)^2-1)*u3_s/6+(norminv(alpha)^3-3*norminv(alpha))*(u4_s-3)/24;
    phi = u1 + sqrt(u2)*phi_s;
    VaR = -(mu+sqrt(htplus)*phi);
end