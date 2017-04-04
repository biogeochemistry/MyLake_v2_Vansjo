% Plots TOC, temperature and oxygen MyLake output vs observations

tlims=[datenum(MyLake_results.m_start):datenum(MyLake_results.m_stop)]';
zinx=find(MyLake_results.z<4.1);

TP_mod = mean((MyLake_results.Pzt(zinx,:)+MyLake_results.PPzt(zinx,:) + MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:)+MyLake_results.DOPzt(zinx,:) +MyLake_results.DOCzt(zinx,:)+MyLake_results.POCzt(zinx,:))', 2);

Chl_mod = mean((MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
Pzt_mod = mean((MyLake_results.Pzt(zinx,:))', 2);
PPzt_mod = mean((MyLake_results.PPzt(zinx,:))', 2);
POCzt_mod = mean((MyLake_results.POCzt(zinx,:))', 2);

load 'obs/store_obs/TOTP.dat' % these are just C&P of vanem ...
load 'obs/store_obs/Cha_aquaM_march_2017.dat' % these are just C&P of vanem ...
load 'obs/store_obs/PO4.dat' % these are just C&P of vanem ...
load 'obs/store_obs/Part.dat' % these are just C&P of vanem ...


startDate = datenum('01-01-2005');
endDate = datenum('12-31-2011');
xData = linspace(startDate,endDate,15);


figure()


subplot(5,1,1)
% dev1 = (PO4_mod)-PO4_obs;
plot(MyLake_results.days, TP_mod,'Color','k','LineWidth',2)
hold on
set(gca,'fontsize',14)
plot(TOTP(:,1), TOTP(:,2), 'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')

ax = gca;
ax.XTick = xData;
datetick(ax, 'x','mmm-yy','keepticks')

ylabel('mg P / m^3')
legend('TIP + DOP (no Chl)')


subplot(5,1,2)
plot(MyLake_results.days, Chl_mod,'Color','k','LineWidth',2)
hold on
plot(Cha_aquaM_march_2017(:,1), Cha_aquaM_march_2017(:,2),'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
ax = gca;
ax.XTick = xData;
datetick(ax, 'x','mmm-yy','keepticks')
set(gca,'fontsize',14)
ylabel('mg P / m^3')
legend('Chl')

subplot(5,1,3)
plot(MyLake_results.days, Pzt_mod,'Color','k','LineWidth',2)
hold on
plot(PO4(:,1), PO4(:,2),'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
ax = gca;
ax.XTick = xData;
datetick(ax, 'x','mmm-yy','keepticks')
set(gca,'fontsize',14)
ylabel('mg P / m^3')
legend('PO4')


subplot(5,1,4)
plot(MyLake_results.days, PPzt_mod+POCzt_mod,'Color','k','LineWidth',2)
hold on
plot(Part(:,1), Part(:,2),'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
ax = gca;
ax.XTick = xData;
datetick(ax, 'x','mmm-yy','keepticks')
set(gca,'fontsize',14)
ylim([0,40])
ylabel('mg P / m^3')
legend('Part')

subplot(5,1,5)
Var = MyLake_results.Szt;
plot(MyLake_results.days, mean((Var(zinx,:))', 2),'Color','k','LineWidth',2)
ax = gca;
ax.XTick = xData;
datetick(ax, 'x','mmm-yy','keepticks')
set(gca,'fontsize',14)
ylabel('mg / m^3')
legend('Any variable')
