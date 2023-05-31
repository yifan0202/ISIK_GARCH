clc;clear;
%% 1. Implied distribution
daylist = dir('../../NminOP_by_day');
ndays = length(daylist)-2;
PARALIST = [];
para0 = [0,log(0.1),0,log(3)];

tic;
for i = 4:ndays+2
    optmat = csvread(['../../NminOP_by_day/',daylist(i).name],0,1);
    rf = optmat(:,6)/100;
    price = optmat(:,10);
    S0 = optmat(:,5);
    K = optmat(:,4);
    T = optmat(:,3)/365;
    CP = optmat(:,2);
    implv = optmat(:,7);
    price = price+(S0-K.*exp(-rf.*T)).*(1-CP);% 看跌期权转化为看涨期权（平价公式)
   
    optmat = [price,S0,T,K,implv,rf];
    Tlist = unique(optmat(:,3));
    paralist = nan(size(Tlist,1),1+1+4);
    disp(daylist(i).name(1:end-4));
    for t = 1:length(Tlist)
        optmat_t = optmat(optmat(:,3)==Tlist(t),:);
        [para] = filter_CFE(optmat_t, para0);
        para0 = para;
        paralist(t,:) = [str2double(daylist(i).name(1:end-4)),Tlist(t),para];
    end
    PARALIST = [PARALIST;paralist];
end
toc;
save('imledmoments_cfe');
function [para] = filter_CFE(optmat,para0)
    b0 = [0,log(mean(sqrt(optmat(:,5)))),0,log(3)];
    options = optimset('Display', 'off','MaxFunEval',15000);
    [para,fval1] = fminsearch(@(x)get_ssr_cfe(optmat,x),b0,options);
    [para2,fval2] = fminsearch(@(x)get_ssr_cfe(optmat,x),para0,options);
%     if fval1<fval2
%         para = para1;
%     else
%         para = para2;
%     end
    para(2) = exp(para(2));
    %para(4) = exp(para(4))+1;
end

function [ssr] = get_ssr_cfe(optmat,b0)
    Ca = optmat(:,1);
    S = optmat(:,2);
    rf = optmat(:,6);
    K = optmat(:,4);
    T = optmat(:,3);
    mu = b0(1);
    sigma = exp(b0(2));
    k3 = b0(3);
    k4 = b0(4);
    
    d = (log(S./K)+mu)/sigma+sigma;
    delta = (mu-rf.*T+sigma^2/2)/sigma;
    
    Sdelta = S.*exp(delta.*sigma);
    Phid = normcdf(d,0,1);
    phid = normcdf(d,0,1);
    
    C = Sdelta.*Phid-K.*exp(-rf.*T).*normcdf(d-sigma,0,1);
    A3 = 1/6*Sdelta.*sigma.*((2*sigma-d).*phid+sigma^2*Phid);
    A4 = 1/24*Sdelta.*sigma.*((d.^2-1-3*sigma*(d-sigma)).*phid+sigma.^3.*Phid);
    A6 = 1/72*Sdelta.*sigma.*(sigma^5.*Phid+(3-6.*d.^2+5.*sigma*(d-(d-sigma).*(sigma.*d-2)-(d-sigma).^3)).*phid);
    
    Y = Ca-C;
    X = [A3,A4,A6];
    ssr = (Y-X*[k3;k4-3;k3^2])'*(Y-X*[k3;k4-3;k3^2]);
end
