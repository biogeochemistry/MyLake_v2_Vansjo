% Plots TOC, temperature and oxygen MyLake output vs observations 

figure(1)

subplot(2,2,1:2)
dev1 = (TOC_mod_matched)-TOC_obs;
plot(Date_mod, (TOC_mod),'Color','r','LineWidth',2)
hold on
plot(Date_matched_toc, TOC_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', Date_mod(1):120:Date_mod(end));
datetick('x','mm/yy')

hold on 

plot (Date_matched_toc(1:length(dev1)),dev1,'Color','k')
plot([Date_matched_toc(1),Date_matched_toc(end)],[ 0 0], 'k-')

hold on 

subplot(2,2,3)

plot(Date_mod_temp, T_mod,'Color','r','LineWidth',2)
hold on
plot(Date_mod_temp, T_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', Date_matched_toc(1):120:Date_matched_toc(end));
datetick('x','mm/yy')

subplot(2,2,4)

plot(Date_mod_O2, O2_mod_top,'Color','r','LineWidth',2)
hold on
plot(Date_mod_O2, O2_obs_top,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
set(gca,'XTick', Date_matched_toc(1):120:Date_matched_toc(end));
datetick('x','mm/yy')
