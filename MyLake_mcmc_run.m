function [nlogl,SS,Nobs] = MyLake_mcmc_run(Xs,partix,fixix,fixval,Ksed,Klake)

% This function distributes the parameters, WITHOUT 10.^trans and
% calculates a negative log likelihood of the pars. I.e. here the values
% are actual and not on log scale.

Pars(partix) = Xs;
Pars(fixix) = fixval;

for pp=1:32 %for each sediment model parameter
    Ksed{pp,1} = Pars(pp);
end

for pp=33:49 %for each lake parameter
    Klake{pp-32,1} = Pars(pp);
end

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 1; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

%% Calls model calibration function
[TP_obs,TP_mod,chl_obs,chl_mod] = fn_MyL_application(Ksed, Klake, use_INCA, run_INCA);

%% calculate innovation of dev1 for each timestep . two cases :  Ice_on or Ice_off
%% updade O2_mod and T_mod  ... According to Honti 

dev1 = TP_mod(:)-TP_obs(:);
dev1(isnan(dev1))=[]; %removing NaNs
SS(1) = nansum((TP_mod(:)-TP_obs(:)).^2);
Nobs(1) = numel(dev1);

dev2 = chl_mod(:)-chl_obs(:);
dev2(isnan(dev2))=[];
SS(2) = nansum((chl_mod(:)-chl_obs(:)).^2);
Nobs(2) = numel(dev2);

nlogl(1) = normlike([0 Pars(50)],dev1); %These are then additive since they are negative log likelihoods
nlogl(2) = normlike([0 Pars(51)],dev2); %These are then additive since they are negative log likelihoods
