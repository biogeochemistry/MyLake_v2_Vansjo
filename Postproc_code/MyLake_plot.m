% Plots TOC, temperature and oxygen MyLake output vs observations

figure()


subplot(4,1,1)
% dev1 = (PO4_mod)-PO4_obs;
plot(TP_date, TP_mod,'Color','k','LineWidth',2)
hold on
plot(TP_date, TP_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', PO4_date(1):120:PO4_date(end));
datetick('x','mm/yy')
ylabel('Concentration, mg P / m^3')
xlabel('Date')
legend('TP')


subplot(4,1,2)

plot(Chl_date, chl_mod,'Color','k','LineWidth',2)
hold on
plot(Chl_date, chl_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', PO4_date(1):120:PO4_date(end));
datetick('x','mm/yy')
ylabel('Concentration, mg P / m^3')
xlabel('Date')
legend('Chl')

subplot(4,1,3)
plot(Part_date, Part_mod,'Color','k','LineWidth',2)
hold on
plot(Part_date, Part_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
ylim([0,40])
set(gca,'XTick', PO4_date(1):120:PO4_date(end));
datetick('x','mm/yy')
ylabel('Concentration, mg P / m^3')
xlabel('Date')
legend('PP')


subplot(4,1,4)

plot(PO4_date, PO4_mod,'Color','k','LineWidth',2')
hold on
plot(PO4_date, PO4_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', PO4_date(1):120:PO4_date(end));
datetick('x','mm/yy')
ylabel('Concentration, mg P / m^3')
xlabel('Date')
legend('PO_4')
