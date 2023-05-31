function p_value = cct(VaR,ret,alpha)
% Input: VaR, ret, alpha
% Output: pvalue of the conditional coverage test
    T = size(ret,1);  
    EFR = sum((ret<=-VaR),1)/T;
    n1 = T*EFR;
    LR_uc = 2*log((EFR^n1*(1-EFR)^(T-n1))/(alpha^n1*(1-alpha)^(T-n1)));
    
    n00 = 0; n01 =0; n10 = 0; n11 =0;
    for t =2:T
        Ut = (ret(t)<= -VaR(t));
        Utm = (ret(t-1)<= -VaR(t-1));
        if (Ut == 0)&&(Utm == 0)
            n00 = n00+1;
        elseif (Ut == 0)&&(Utm ==1)
            n10 = n10+1;
        elseif (Ut == 1)&&(Utm == 0)
            n01 = n01+1;
        else
            n11 = n11+1;
        end
    end
    
    pi01 = n01/(n00+n01);
    pi11 = n11/(n10+n11);
    
    LR_cc = LR_uc+2*log((pi01^n01*(1-pi01)^n00*pi11^n11*(1-pi11)^n10)/(EFR^(n01+n11)*(1-EFR)^(n00+n10)));
    
    p_value = 1-chi2cdf(LR_cc,2);
end