%% plots the diagnostic graphs for residuals of TOC, and O2 (with and without ice)
%% carbon
figure('name','Residual diagnostic for TOC'); 

subplot(2,2,1) % time series of obs, sim and residuals
    plot(Date_matched_toc, (TOC_mod_matched),'Color','r','LineWidth',2)
        hold on
    plot(Date_matched_toc, TOC_obs,'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
    set(gca,'XTick', Date_mod(1):120:Date_mod(end));
    datetick('x','mm/yy')
        hold on
    dev1 = (TOC_mod_matched)-TOC_obs;
    plot (Date_matched_toc(1:length(dev1)),dev1,'Color','k')
    plot([Date_matched_toc(1),Date_matched_toc(end)],[ 0 0], 'k-')
        hold on
subplot(2,2,2) % residuals vs predicted values
    temp = TOC_obs ;
    temp(isnan(TOC_obs))=[]; % removing nans from TOC_obs
    dev1(isnan(dev1))=[]; %removing NaNs
    plot(dev1, temp,'Color','r','LineWidth',2)
        hold on
subplot(2,2,3) % autocorrelation of residuals
    autocorr(dev1)
        hold on 
subplot(2,2,4) % QQ plots 
    qqplot(dev1)

%% oxygen ice_free period
figure('name','Residual diagnostic for DO ice-free period'); 

subplot(2,2,1) % time series of obs, sim and residuals
    plot(Date_mod_O2(Ice_free_vec), O2_mod(Ice_free_vec),'Color','r','LineWidth',2)
        hold on
    plot(Date_mod_O2(Ice_free_vec), O2_obs(Ice_free_vec),'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
    set(gca,'XTick', Date_matched_toc(1):120:Date_matched_toc(end));
    datetick('x','mm/yy')
        hold on 
    plot (Date_mod_O2(Ice_free_vec),dev2,'Color','k')
    plot([Date_mod_O2(1),Date_mod_O2(end)],[ 0 0], 'k-')
    hold on 

subplot(2,2,2) % residuals vs predicted values
    temp = O2_obs(Ice_free_vec) ;
    temp(isnan(O2_obs(Ice_free_vec)))=[]; % removing nans from O2_obs
    nan_dev2 = dev2;
    nan_dev2(isnan(dev2))=[]; %removing NaNs
    plot(nan_dev2, temp,'Color','r','LineWidth',2)
        hold on

subplot(2,2,3) % autocorrelation of residuals
autocorr(nan_dev2)
    hold on 

subplot(2,2,4) % QQ plots 
qqplot(nan_dev2)

%% oxygen covered period
figure('name','Residual diagnostic for DO ice-covered period'); 

subplot(2,2,1) % time series of obs, sim and residuals
    plot(Date_mod_O2(Ice_cover_vec), O2_mod(Ice_cover_vec),'Color','r','LineWidth',2)
        hold on
    plot(Date_mod_O2(Ice_cover_vec), O2_obs(Ice_cover_vec),'Color','b','marker','diamond','MarkerSize',4,'LineStyle','none')
    set(gca,'XTick', Date_matched_toc(1):120:Date_matched_toc(end));
    datetick('x','mm/yy')
        hold on 
    plot (Date_mod_O2(Ice_cover_vec),dev3,'Color','k')
    plot([Date_mod_O2(1),Date_mod_O2(end)],[ 0 0], 'k-')
    hold on 

subplot(2,2,2) % residuals vs predicted values
    temp = O2_obs(Ice_cover_vec) ;
    temp(isnan(O2_obs(Ice_cover_vec)))=[]; % removing nans from O2_obs
    nan_dev3 = dev3;
    nan_dev3(isnan(dev3))=[]; %removing NaNs
    plot(nan_dev3, temp,'Color','r','LineWidth',2)
        hold on

subplot(2,2,3) % autocorrelation of residuals
autocorr(nan_dev3)
    hold on 

subplot(2,2,4) % QQ plots 
qqplot(nan_dev3)