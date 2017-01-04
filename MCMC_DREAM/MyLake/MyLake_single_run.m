clear all; 					% range
tic
MyLake_setPriors_for_MCMC % prepares parameter files

m_start=[1984, 5, 31]; % 1973 1 1 is the earliest possible given the provided input files8
m_stop=[1986, 12, 31]; % 2012, 12, 31 is the latest possible given the provided input files  

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

[TP_obs,TP_mod,chl_obs,chl_mod] = fn_MyL_application(m_start, m_stop, K_values_sediments, K_values_lake, use_INCA, run_INCA); % runs the model and outputs obs and sim

dev1 = TP_mod(:)-TP_obs(:);
dev1(isnan(dev1))=[]; %removing NaNs
SS(1) = nansum((TP_mod(:)-TP_obs(:)).^2);
Nobs(1) = numel(dev1);

dev2 = chl_mod(:)-chl_obs(:);
dev2(isnan(dev2))=[];
SS(2) = nansum((chl_mod(:)-chl_obs(:)).^2);
Nobs(2) = numel(dev2);

clear dev1 dev2 Nobs
toc