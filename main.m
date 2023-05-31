clc;clear;
%% 0.1. Load Data
[Price,Date] = xlsread('../Data/50ETF.xlsx');
Date = Date(6:end,1);
Date = datestr(Date,'yyyymmdd'); 
Date = str2num(Date);
ret = log(Price(2:end))-log(Price(1:end-1));
Date = Date(2:end);
startdate = 20150210;
enddate = 20181228;
ret = ret((Date>=startdate)&(Date<=enddate),:);
DateUse = Date((Date>=startdate)&(Date<=enddate),:);

load(['1.Implied_moments/data.mat']);
ind = ismember(DateUse,outmat30(:,1));
miss = DateUse(ind==0,:);
DateUse = DateUse(ind == 1,:);
ret = ret(ind == 1,:);

clear Date Price ind
%% 0.2. Summary Stat
data = [ret*100, outmat30(:,4:5)];
stat = [mean(data)', median(data)', std(data)', skewness(data)',kurtosis(data)'-3, min(data)', max(data)'];
xlswrite('../Data/Tables.xlsx',stat,'summary_stat','c3');
%% 1.1. GARCH
cd('2.GARCH')
options = optimset('Display', 'off','MaxFunEval',15000);
[GJR.para,~,GJR.ht, GJR.VCVrobust] = tarch(ret,1,1,1,'NORMAL',[],[],options);
GJR.tstat = GJR.para./sqrt(diag(GJR.VCVrobust));%杠杆效应不显著转为用GARCH（1，1）
% [G.para,G.LL,G.ht, G.VCVrobust] = tarch(ret,1,0,1,'NORMAL',[],[],options);
% G.tstat = G.para./sqrt(diag(G.VCVrobust));
[GH.para, GH.ht, GH.LL, GH.se, GH.htplus] = GARCH(ret);
GH.tstat = GH.para'./GH.se;
cd('..');
%% 1.2. GARCH-SK
cd('3.GARCH_SK')
[SK.para, SK.ht, SK.zt,SK.s,SK.k,SK.LL, SK.LL_r,SK.se, SK.htplus,SK.splus,SK.kplus] = GARCH_SK(ret);
SK.tstat = SK.para'./SK.se;
cd('..');
%% 1.3. GARCH-ISIK
cd('4.GARCH_ISIK')
[ISIK.para, ISIK.ht,ISIK.zt,ISIK.s,ISIK.k,ISIK.eta,ISIK.vartheta, ISIK.LL,ISIK.LL_r,ISIK.LL_IS,ISIK.LL_IK, ISIK.se, ISIK.htplus,ISIK.splus,ISIK.kplus] = GARCH_ISIK(ret, outmat30(:,4:5));
ISIK.tstat = ISIK.para'./ISIK.se;

fprintf('hat_beta2:%f\n',ISIK.para(7)*ISIK.para(10));
fprintf('delta2:%f\n',SK.para(7));
fprintf('\nhat_gamma2:%f\n',ISIK.para(13)*ISIK.para(16));
fprintf('theta2:%f\n',SK.para(10));
cd('..');
%% 1.4. Table
outmat1 = nan(16,2);
outmat2 = nan(16,2);
outmat3 = nan(16,2);
outmat1(1:4,1) = GH.para';
outmat1(1:4,2) = GH.se;
outmat2(:,1) = [SK.para(1,1:6)';nan(3,1);SK.para(1,7:9)';nan(3,1);SK.para(1,10)'];
outmat2(:,2) = [SK.se(1:6,1);nan(3,1);SK.se(7:9,1);nan(3,1);SK.se(10,1)];
outmat3(:,1) = ISIK.para(1,1:16)';
outmat3(:,2) = ISIK.se(1:16,1);
outmat1 = reshape(outmat1',32,1);
outmat2 = reshape(outmat2',32,1);
outmat3 = reshape(outmat3',32,1);

outmat1(33) = GH.LL;
outmat2(33) = SK.LL_r;
outmat3(33) = ISIK.LL_r;
outmat2(34) = SK.LL;
outmat3(34) = ISIK.LL;
outmat1 = [outmat1(1:8);nan;outmat1(9:20);nan;outmat1(21:33)];
outmat2 = [outmat2(1:8);nan;outmat2(9:20);nan;outmat2(21:34)];
outmat3 = [outmat3(1:8);nan;outmat3(9:20);nan;outmat3(21:34)];

output = '../Data/Tables.xlsx';
output_sheet = 'insample';
xlswrite(output,outmat1,output_sheet,'c4');
xlswrite(output,outmat2,output_sheet,'d4');
xlswrite(output,outmat3,output_sheet,'e4');
%% 2. Correct moments
[ISIK.u1, ISIK.u2, ISIK.u3, ISIK.u4] = Transf(ISIK.s,ISIK.k);
TT.u3 = ISIK.u3;
TT.u4 = ISIK.u4;
TT.s = ISIK.s;
TT.k = ISIK.k;

figure(1)
subplot(2,1,1)
plot(TT.Time,TT.u3);
hold on 
plot(TT.Time,TT.s);
legend({'u_{3t}','s_t'});
subplot(2,1,2)
plot(TT.Time,TT.u4);
hold on
plot(TT.Time,TT.k);
legend({'u_{4t}','k_t'});

save data
%% 3.1. OOS forecast
% cd('5.GARCH_t')
% [GHt.para, GHt.ht, GHt.LL, GHt.se, GHt.htplus] = GARCH_t(ret);
% cd('..');
% [Gt.para,Gt.LL,Gt.ht, Gt.VCVrobust] = tarch(ret,1,0,1,'STUDENTST',[],[],options);
tic;
[VaR_ISIK,VaR_SK,VaR_G,VaR_Gt] = VaR_OOS(ret, outmat30(:,4:5),250);
toc;
%% 3.2. EFR & the unconditional coverage tests & conditional coverage tests
[Test_ISIK,Test_SK,Test_G,Test_Gt] = VaRbt(VaR_ISIK.VaR,VaR_SK.VaR,VaR_G.VaR,VaR_Gt.VaR,ret,250);

output = '../Data/Tables.xlsx';
output_sheet = 'Test';
xlswrite(output,100*Test_ISIK.EFR',output_sheet,'d5');
xlswrite(output,100*Test_SK.EFR',output_sheet,'e5');
xlswrite(output,100*Test_G.EFR',output_sheet,'f5');
xlswrite(output,100*Test_Gt.EFR',output_sheet,'g5');

xlswrite(output,Test_ISIK.p_uc',output_sheet,'d12');
xlswrite(output,Test_SK.p_uc',output_sheet,'e12');
xlswrite(output,Test_G.p_uc',output_sheet,'f12');
xlswrite(output,Test_Gt.p_uc',output_sheet,'g12');

xlswrite(output,Test_ISIK.p_cc',output_sheet,'d19');
xlswrite(output,Test_SK.p_cc',output_sheet,'e19');
xlswrite(output,Test_G.p_cc',output_sheet,'f19');
xlswrite(output,Test_Gt.p_cc',output_sheet,'g19');

[Test_ISIK_sk,Test_SK_sk,~,~] = VaRbt(VaR_ISIK.VaR_sk,VaR_SK.VaR_sk,VaR_G.VaR,VaR_Gt.VaR,ret,250);
output = '../Data/Tables.xlsx';
output_sheet = 'Test_sk';
xlswrite(output,100*Test_ISIK_sk.EFR,output_sheet,'d5');
xlswrite(output,100*Test_SK_sk.EFR,output_sheet,'d6');
xlswrite(output,Test_ISIK_sk.p_uc,output_sheet,'d8');
xlswrite(output,Test_SK_sk.p_uc,output_sheet,'d9');
%% 4.Plot
TT.ret = ret;
for i =1:6
    eval(['TT.VaR_G',num2str(i),'=-VaR_G.VaR(:,',num2str(i),');']);
    eval(['TT.VaR_Gt',num2str(i),'=-VaR_Gt.VaR(:,',num2str(i),');']);
    eval(['TT.VaR_ISIK',num2str(i),'=-VaR_ISIK.VaR(:,',num2str(i),');']);
    eval(['TT.VaR_SK',num2str(i),'=-VaR_SK.VaR(:,',num2str(i),');']);
end
T = TT(~isnan(TT.VaR_G1),:);

figure(2)
subplot(3,2,1)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G1);
hold on
plot(T.Time,T.VaR_Gt1);
hold on
plot(T.Time,T.VaR_ISIK1);
hold on
plot(T.Time,T.VaR_SK1);
h1 = legend({'Return','VaR(0.5%):GARCH','VaR(0.5%):GARCH-t','VaR(0.5%):GARCH-ISIK','VaR(0.5%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

subplot(3,2,2)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G2);
hold on
plot(T.Time,T.VaR_Gt2);
hold on
plot(T.Time,T.VaR_ISIK2);
hold on
plot(T.Time,T.VaR_SK2);
h1 = legend({'Return','VaR(1%):GARCH','VaR(1%):GARCH-t','VaR(1%):GARCH-ISIK','VaR(1%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

subplot(3,2,3)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G3);
hold on
plot(T.Time,T.VaR_Gt3);
hold on
plot(T.Time,T.VaR_ISIK3);
hold on
plot(T.Time,T.VaR_SK3);
h1 = legend({'Return','VaR(1.5%):GARCH','VaR(1.5%):GARCH-t','VaR(1.5%):GARCH-ISIK','VaR(1.5%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

subplot(3,2,4)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G4);
hold on
plot(T.Time,T.VaR_Gt4);
hold on
plot(T.Time,T.VaR_ISIK4);
hold on
plot(T.Time,T.VaR_SK4);
h1 = legend({'Return','VaR(2%):GARCH','VaR(2%):GARCH-t','VaR(2%):GARCH-ISIK','VaR(2%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

subplot(3,2,5)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G5);
hold on
plot(T.Time,T.VaR_Gt5);
hold on
plot(T.Time,T.VaR_ISIK5);
hold on
plot(T.Time,T.VaR_SK5);
h1 = legend({'Return','VaR(2.5%):GARCH','VaR(2.5%):GARCH-t','VaR(2.5%):GARCH-ISIK','VaR(2.5%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

subplot(3,2,6)
plot(T.Time,T.ret);
hold on 
plot(T.Time,T.VaR_G6);
hold on
plot(T.Time,T.VaR_Gt6);
hold on
plot(T.Time,T.VaR_ISIK6);
hold on
plot(T.Time,T.VaR_SK6);
h1 = legend({'Return','VaR(5%):GARCH','VaR(5%):GARCH-t','VaR(5%):GARCH-ISIK','VaR(5%):GARCH-SK'},'Location','SouthOutside');
set(h1,'Orientation','horizon')
h1.ItemTokenSize = [8,8];

savefig();