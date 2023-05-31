function [ISIK,SK,G,Gt] = VaRbt(VaR_ISIK,VaR_SK,VaR_G,VaR_Gt,ret,windowsize)
% Input: VaR_ISIK,VaR_SK,VaR_G,VaR_Gt
%           ret(return data)
% Ouput: ISIK, SK, G, Gt( structure, including EFR, p-values for the
% unconditional coverage tests and conditional coverage tests
    VaR_ISIK = VaR_ISIK(~isnan(VaR_ISIK(:,1)),:);
    VaR_SK = VaR_SK(~isnan(VaR_SK(:,1)),:);
    VaR_G = VaR_G(~isnan(VaR_G(:,1)),:);
    VaR_Gt = VaR_Gt(~isnan(VaR_Gt(:,1)),:);
    ret = ret(windowsize+1:end);
    n = size(VaR_G,2);
    T = size(ret,1);
    
    EFR_ISIK = sum((repmat(ret,1,n)<=-VaR_ISIK),1)/T;
    EFR_SK = sum((repmat(ret,1,n)<=-VaR_SK),1)/T;
    EFR_G = sum((repmat(ret,1,n)<=-VaR_G),1)/T;
    EFR_Gt = sum((repmat(ret,1,n)<=-VaR_Gt),1)/T;
    
    alpha = [0.5/100,1/100,1.5/100,2/100,2.5/100,5/100];
    
    uct_ISIK = nan(1,n);
    uct_SK = nan(1,n);
    uct_G = nan(1,n);
    uct_Gt = nan(1,n);
    
    cct_ISIK = nan(1,n);
    cct_SK = nan(1,n);
    cct_G = nan(1,n);
    cct_Gt = nan(1,n);
    
    for i =1:n
        uct_ISIK(i) = uct(EFR_ISIK(i),alpha(i),T);
        uct_SK(i) = uct(EFR_SK(i),alpha(i),T);
        uct_G(i) = uct(EFR_G(i),alpha(i),T);
        uct_Gt(i) = uct(EFR_Gt(i),alpha(i),T);
        
        cct_ISIK(i) = cct(VaR_ISIK(:,i),ret,alpha(i));
        cct_SK(i) = cct(VaR_SK(:,i),ret,alpha(i));
        cct_G(i) = cct(VaR_G(:,i),ret,alpha(i));
        cct_Gt(i) = cct(VaR_Gt(:,i),ret,alpha(i));
    end
    ISIK.EFR = EFR_ISIK;
    SK.EFR = EFR_SK;
    G.EFR = EFR_G;
    Gt.EFR = EFR_Gt;
    
    ISIK.p_uc = uct_ISIK;
    SK.p_uc = uct_SK;
    G.p_uc = uct_G;
    Gt.p_uc = uct_Gt;
    
    ISIK.p_cc = cct_ISIK;
    SK.p_cc = cct_SK;
    G.p_cc = cct_G;
    Gt.p_cc = cct_Gt;
end

