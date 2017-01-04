% Plots oxygen and temperature at 2 depths
% Start_date = datenum(2010,5,27);
% Date_ticks = [1:1315] + Start_date;
% Date_string = datestr(Date_ticks, 'dd.mm.yyyy');
%
figure

subplot(2,1,1)
plot(Date_mod_O2,O2_mod_top,'--r','LineWidth',1.5)
hold on
plot(Date_mod_O2,O2_obs_top,'LineWidth',1.5)
title('Oxygen modeled (red) and observed (blue) at 1 m depth')
ylabel('O2 concentration (mg/m3)')
% xlabel('Simulation date from 2010/04 to 2012-06 (mm)')
datetick('x','mmm')

subplot(2,1,2)
plot(Date_mod_O2,O2_mod_bottom,'--r','LineWidth',1.5)
hold on
plot(Date_mod_O2,O2_obs_bottom,'LineWidth',1.5)
title('Oxygen modeled (red) and observed (blue) at 9 m depth')
ylabel('O2 concentration (mg/m3)')
% xlabel('Simulation date from 2010/04 to 2012-06 (month)')
datetick('x','mmm')

%
% figure

% subplot(2,1,1)
% plot(Date_mod_O2,T_mod,'--r','LineWidth',1.5)
% hold on
% plot(Date_mod_temp,T_obs,'LineWidth',1.5)
% title('Temperature modeled (red) and observed (blue) at 1 m depth')
% ylabel('Temperature (C)')
% % xlabel('Simulation date from 2010/04 to 2012-06 (mm)')
% datetick('x','mmm')

% subplot(2,1,2)
% plot(Date_mod_temp,T_mod,'--r','LineWidth',1.5)
% hold on
% plot(Date_mod_O2,T_obs,'LineWidth',1.5)
% title('Temperature modeled (red) and observed (blue) at 9 m depth')
% ylabel('Temperature (C)')
% % xlabel('Simulation date from 2010/04 to 2012-06 (mm)')
% datetick('x','mmm')
