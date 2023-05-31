function [ISIK,SK,GH,GHt] = VaR_OOS(ret, SK_data, windowsize)
% Input: ret(return); ISIK(implied skewness and kurtosis)
% Output: VaR(0.5, 1.0, 1.5, 2.0, 2.5, 5.0%)
    T = size(ret,1);
    VaR_ISIK = nan(T,6);  
    VaR_ISIK_sk = nan(T,6);
    VaR_SK = nan(T,6);
    VaR_SK_sk = nan(T,6);
    VaR_GH = nan(T,6);   
    VaR_GHt = nan(T,6);
    
    ll_ISIK = nan(T,1);
    ll_SK = nan(T,1);
    ll_GH = nan(T,1);
    ll_GHt = nan(T,1);
    
    para_ISIK = nan(T,20);
    para_SK = nan(T,12);
    para_GH = nan(T,4);
    para_GHt = nan(T,5);
   
    h_ISIK = nan(T,1);
    h_SK = nan(T,1);
    h_GH = nan(T,1);
    h_GHt = nan(T,1);
    
    alpha = [0.5/100, 1/100, 1.5/100, 2/100, 2.5/100, 5/100];
    
    for t = windowsize+1:T
        insample = ret(t-windowsize:t-1);
        ISIK_mat = SK_data(t-windowsize:t-1,:);
        cd('4.GARCH_ISIK')
        [ISIK_para,~,~,~,~,~,~,ISIK_ll,~,~,~,~, ISIK_htplus,ISIK_splus,ISIK_kplus] = GARCH_ISIK(insample, ISIK_mat);
        cd('..');
        cd('3.GARCH_SK')
        [SK_para,~,~,~,~,SK_ll,~,~,SK_htplus,SK_splus,SK_kplus] = GARCH_SK(insample);
        cd('..');
        [ISIK_u1, ISIK_u2, ISIK_u3, ISIK_u4] = Transf(ISIK_splus,ISIK_kplus);
        ISIK_u = [ISIK_u1, ISIK_u2, ISIK_u3, ISIK_u4];
        [SK_u1, SK_u2, SK_u3, SK_u4] = Transf(SK_splus,SK_kplus);
        SK_u = [SK_u1, SK_u2, SK_u3, SK_u4];    
        cd('2.GARCH')
        [GH_para,~,GH_ll,~,GH_htplus] = GARCH(insample);
        cd('..');
        cd('5.GARCH_t')
        [GHt_para,~,GHt_ll,~,GHt_htplus] = GARCH_t(insample);
        cd('..');  
        for i = 1:size(alpha,2)
            VaR_ISIK(t,i) = getVaR(ISIK_para(1),ISIK_u, alpha(i),ISIK_htplus);
            VaR_SK(t,i) = getVaR(SK_para(1),SK_u, alpha(i),SK_htplus);
            VaR_ISIK_sk(t,i) = getVaR_sk(ISIK_para(1),ISIK_splus,ISIK_kplus,alpha(i),ISIK_htplus);
            VaR_SK_sk(t,i) = getVaR_sk(SK_para(1),SK_splus,SK_kplus,alpha(i),SK_htplus);
            VaR_GH(t,i) = -(GH_para(1)+sqrt(GH_htplus)*norminv(alpha(i),0,1));
            VaR_GHt(t,i) = -(GHt_para(1)+sqrt(GHt_htplus)*tinv(alpha(i),GHt_para(end)));
        end
        ll_ISIK(t) = ISIK_ll;
        ll_SK(t) = SK_ll;
        ll_GH(t) = GH_ll;
        ll_GHt(t) = GHt_ll;
        
        para_ISIK(t,:) = ISIK_para;
        para_SK(t,:) = SK_para;
        para_GH(t,:) = GH_para;
        para_GHt(t,:) = GHt_para;
        
        h_ISIK(t) = ISIK_htplus;
        h_SK(t) = SK_htplus;
        h_GH(t) = GH_htplus;
        h_GHt(t) = GHt_htplus;
        fprintf('%i\n',t);
    end
    ISIK.VaR = VaR_ISIK;
    ISIK.ll =ll_ISIK;
    ISIK.para = para_ISIK;
    ISIK.VaR_sk = VaR_ISIK_sk;
    ISIK.htplus = h_ISIK;
    
    SK.VaR = VaR_SK;
    SK.ll =ll_SK;
    SK.para = para_SK; 
    SK.VaR_sk = VaR_SK_sk;
    SK.htplus = h_SK;
    
    GH.VaR = VaR_GH;
    GH.ll =ll_GH;
    GH.para = para_GH;
    GH.htplus = h_GH;
    
    GHt.VaR = VaR_GHt;
    GHt.ll =ll_GHt;
    GHt.para = para_GHt;
    GHt.htplus = h_GHt;
end