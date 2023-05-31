clc;clear;

load('imledmoments_cfe');
daylist = unique(PARALIST(:,1));
T = length(daylist);
outmat30 = nan(T,1+4);

for t = 1:T
    idx = (PARALIST(:,1) == daylist(t));
    outmat30(t,1) = daylist(t);
    if sum(idx)>2
        daysec = PARALIST(idx,:);
        daysec = sortrows(daysec,2);
        if daysec(1,2)>30/365
            outmat30(t,2:5) = daysec(1,3:6);
        elseif daysec(end,2)<30/365
            outmat30(t,2:5) = daysec(end,3:6);
        else
            for i = 2:5
                outmat30(t,i) = pchip(daysec(:,2),daysec(:,i+1),30/365);
            end
         end
    end
end

date=datetime(outmat30(:,1),'ConvertFrom','yyyymmdd');
TT = timetable(outmat30(:,2),outmat30(:,3),outmat30(:,4),outmat30(:,5),'RowTimes',date,'VariableNames',{'mu','sigma','k3','k4'});
%转化为日度
% TT.mu = TT.mu/30;
% TT.sigma = TT.sigma/30;
% TT.k3 = TT.k3/30;
% TT.k4 = TT.k4/30;

figure(1)
subplot(2,2,1)
plot(TT.Time,TT.mu);
legend({'\mu'});
subplot(2,2,2)
plot(TT.Time,TT.sigma);
legend({'\sigma'});
subplot(2,2,3)
plot(TT.Time,TT.k3);
legend({'Skewness'});
subplot(2,2,4)
plot(TT.Time,TT.k4);
legend({'Kurtosis'});

figure(2)
subplot(2,1,1)
plot(TT.Time,TT.k3);
hold on
plot(TT.Time, zeros(length(TT.Time)),'--','Color',[0.85,0.33,0.10]);
legend({'Implied Skewness','Skewness of Normal Distribution'});
subplot(2,1,2);
plot(TT.Time,TT.k4);
hold on
plot(TT.Time, 3*ones(length(TT.Time)),'--','Color',[0.85,0.33,0.10]);
legend({'Implied Skewness','Kurtosis of Normal Distribution'});

save data outmat30 TT
