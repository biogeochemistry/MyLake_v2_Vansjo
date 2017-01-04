load('sediment_data_basin2.mat');
load('MyLake_results_basin2.mat');

Fe2 = [MyLake_results_basin2{19, 1}/55845; sediment_data_basin2{5, 1}];

% NO3 = [MyLake_results_basin2{14, 1}/62004; sediment_data_basin2{19, 1}];

% NH4 = [MyLake_results_basin2{15, 1}/18038; sediment_data_basin2{20, 1}];

O2 = [MyLake_results_basin2{13, 1}/31998; sediment_data_basin2{1, 1}];

O2_flux_WC = MyLake_results_basin2{15, 1};

O2_flux_sed = sediment_data_basin2{42, 1}{1, 1};

pH_sed = sediment_data_basin2{31, 1};

start_date = datenum(sediment_data_basin2{39,1}{2,2});


% f1 = figure(1)
% mesh(O2)
% % h_legend=legend('$NO_3$ concentration');
% % set(h_legend,'Interpreter','latex')
% set(f1, 'defaulttextinterpreter','latex');
% % tit = title('$ NO_3 $ concentration with bioirrigation');
% xlabel('time, $t$ $ \left[ day \right] $');
% ylabel('depth mesh point, $n$');
% zlabel('$[C]$ , $\frac{\mu mol}{cm^3}$');
% % set(h_legend,'FontSize',24);
% set(gca,'FontSize',24)
% set(gca, 'XMinorTick', 'on')
% set(gca, 'YMinorTick', 'on')
% set( findobj(gca,'type','line'), 'LineWidth', 2);



f1 = figure(2)
plot((1:size(O2_flux_WC,2))+start_date,O2_flux_WC/32, (1:size(O2_flux_sed,2)) + start_date, -O2_flux_sed/32)
h_legend=legend('by WC','by Sediments');
% set(h_legend,'Interpreter','latex')
set(f1, 'defaulttextinterpreter','latex');
% tit = title('$ NO_3 $ concentration with bioirrigation');
xlabel('time, $t$ $ \left[ day \right] $');
ylabel('$[O_2]$ production, $\frac{mmol}{m^2 \cdot day}$' );
% zlabel('$[C]$ , $\frac{\mu mol}{cm^3}$');
% set(h_legend,'FontSize',24);
set(gca,'FontSize',24)
set(gca, 'XMinorTick', 'on')
set(gca, 'YMinorTick', 'on')
set( findobj(gca,'type','line'), 'LineWidth', 2);
datetick('x','mm/yy')

% f1 = figure(3)
% mesh(pH_sed)
% % h_legend=legend('$NO_3$ concentration');
% % set(h_legend,'Interpreter','latex')
% set(f1, 'defaulttextinterpreter','latex');
% % tit = title('$ NO_3 $ concentration with bioirrigation');
% xlabel('time, $t$ $ \left[ day \right] $');
% ylabel('depth mesh point, $n$');
% zlabel('pH');
% % set(h_legend,'FontSize',24);
% set(gca,'FontSize',24)
% set(gca, 'XMinorTick', 'on')
% set(gca, 'YMinorTick', 'on')
% set( findobj(gca,'type','line'), 'LineWidth', 2);