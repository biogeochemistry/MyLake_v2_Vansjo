%% Initialisation 
clear all

cd ../MyLake 

addpath(genpath('C:\OneDrive\PROJECT\0-12111_MARS\4.4_Northern_basin\Vansjø'))

tic
total_runs = 0;
%% Generatiing Priors for MyLake 

MyLake_setPriors_for_MCMC
% Setting a structrued array with priors, Not for the sigmas.

%& Main model call

% [T_obs,T_mod,O2_obs,O2_mod] = Calibration(K_values_sediments, K_values_lake);

% where K_values_sediments are parameters (1:32) and K_values_lake are
% parameters 33:49
MCMCPar.n           = 8;       % Number of chains
MCMCPar.b           = 0.1;      % For drawing e
MCMCPar.b_star      = 1e-3;     % For drawing eps
MCMCPar.CRupdate    = 'Yes';    % Updating the crossover values or not
MCMCPar.nCR         =  3;       % Number of crossover values ... 
MCMCPar.delta       =  1;       % Number of different pairs of chains to look at. set to 0 original 3Must be equal or less than (MCMCPar.n-1)/2
MCMCPar.gammaisone  =  10;      % Setting gamma==1 for every t=5
MCMCPar.lhsInit     = 'Yes';    % Initialize by latin hypercude sampling
MCMCPar.checkGR     =   10;     % Checking Gelman-Rubin diagnostic every t. Not before burnin is over
MCMCPar.checkoutlier=   25;     % Checking for outlier every t=. Not after burnin is over. Not implemented for the testing here, since likelihoods are not calculated.
MCMCPar.tmax        = 500;        % Stop after this many iterations
MCMCPar.burnin      = 100;     % End of burninperiod.
plt                 = 100;       %Plot every plth iteration
outlierc 			= 0;		% To count the number of outliers removed during burnin. 

Priors(50).name = {'Error variance TP'};
Priors(51).name = {'Error variance Chl'};
Priors(50).min  = 1e-3;
Priors(50).max  = 5e3;
Priors(51).min  = 1e-3;
Priors(51).max  = 1e2;
%%
% Indexing which parameters to vary. I there is no interval, then no
% variation. So clever ! 
tix = 1;
tix2= 1;

for ii=1:51% should be generalized ... lenght of whatever % 51 because its accounting for the 2 SS
    if (Priors(ii).max-Priors(ii).min)>0
        Priors(ii).xinx = tix;
        partix(tix)=ii;
        MinMax(tix,:) = [Priors(ii).min Priors(ii).max];
        tix = tix+1;
    else
        fixix(tix2) = ii;
        fixval(tix2) = Priors(ii).min;
        tix2=tix2+1;
    end
end
MCMCPar.nparam = numel(partix);
MinMax = log10(MinMax);
%%
% Some initial testing;
% figure(1)
% subplot(2,1,1),hist(O2_mod'-O2_obs')
% subplot(2,1,2),hist(T_mod'-T_obs')
% % For initial try assuming that we have 2 error variances, one for temp and
% % one for o2
% dev1 = O2_mod(:)-O2_obs(:);
% dev1(isnan(dev1))=[]; %removing NaNs
% dev2 = T_mod(:)-T_obs(:);
% dev2(isnan(dev2))=[];
%
% nlogL = normlike([0 1e-2],dev2) %These are then additive since they are negative log likelihoods


%% Priorthoughts
% I think we might set all the priors to be loguniform, since they are rate
% constants mostly and have ranges of type 1e-4 1e3. We can thus treat all
% as uniform on the log scale, also the error variances.
% We will thus sample all parameters on the log scale
%% Initializing
%=====================================================
%Initializing the CR-values
if strcmp(MCMCPar.CRupdate,'Yes') %if updating the cr values
    L = ones(1,MCMCPar.nCR); %L values used %[250811,JOS] now using ones since it will become Inf if /0
    CRv = (1:MCMCPar.nCR)./MCMCPar.nCR;  %initial CR values.
    pCR = ones(1,MCMCPar.nCR).*(1./MCMCPar.nCR); %prob for initial CR values
    deltas = zeros(MCMCPar.nCR,1); %for tracking individual delta values
    delta_tot = 0;
end

Chains = NaN*ones(MCMCPar.nparam,MCMCPar.n,MCMCPar.tmax); %chains now.
Accept = NaN*ones(MCMCPar.tmax,1);
nLogHist = NaN*ones(MCMCPar.n,MCMCPar.tmax);
SS_Hist = NaN*ones(MCMCPar.n,2,MCMCPar.tmax); %assuming two SS for each chain per iter.
GR_Hist = NaN*ones(MCMCPar.nparam,ceil(MCMCPar.tmax/MCMCPar.checkGR));
gr_tix = 1;
Accept(1) = 1;
%%
%Init = rand(MCMCPar.n,MCMCPar.nparam)';
Init = lhsdesign(MCMCPar.n,MCMCPar.nparam)';
% TRANSPOSING THIS TO MAKE IT (params, chains)
% All are LOG scale
for cc=1:MCMCPar.n %for each chain
    % Xs(:,cc)  = 10.^(log10(MinMax(:,1)) + (Init(cc,:)').*(log10(MinMax(:,2))-log10(MinMax(:,1))))
    Xs(:,cc)  = ((MinMax(:,1)) + (Init(:,cc)).*((MinMax(:,2))-(MinMax(:,1))));
    % THESE Xs are on a log scale, translate with 10.^ to generate parameter
    % values actually used.
end
% X(42 and 43) are the sigmas.

 X_actual = 10.^(Xs);
% parfor sometimes crashes, sometimes not. mystery.  
    disp('Initializing chains ... ')
        
for cc=1:MCMCPar.n %for all chains, first run.
   [nlogl(cc,:),SS(cc,:),Nobs] = MyLake_mcmc_run(X_actual(:,cc),partix,fixix,fixval,K_values_sediments,K_values_lake)
   %disp('cooling down')
   %pause(2)
end


%% Parallel model runs and iterations
nLogHist(:,1) = sum(nlogl,2); %adding them up
SS_Hist(:,:,1) = SS;
X_now = Xs;
Chains(:,:,1) = Xs;
disp('Starting iterations ... ')
tt=1;
while tt< MCMCPar.tmax %Main iteration loop
    
    tt=tt+1;
    % Put in proposals here
    [L,X_prop,WhichM] = generateProposalv4(MCMCPar,X_now,L,tt,pCR);
    % Now sigmas are as the other ones, i.e. MCMC updateed and not Gibbs.
    % This can be changed, but then all arrays will need to change.
    
    %
    % Evaluate priors
    ignorenow = []; %for collecting which runs to ignore
    SS_new = NaN*ones(MCMCPar.n,2);
    nlogl_new = NaN*ones(MCMCPar.n,2);
    alphas = NaN*ones(MCMCPar.n,1);
    % OLD testing for outside parameter bounds. Taken out JOS060514
    %     for cc=1:MCMCPar.n
    %         if any(X_prop(:,cc)<MinMax(:,1))>0 | any(X_prop(:,cc)>MinMax(:,2))>0 %if any are outside
    %             SS_new(cc,:) = Inf;
    %             nlogl_new(cc,:)    = NaN; %???
    %             ignorenow = [ignorenow;cc];
    %         end
    %     end
    
    % NEW (JOS060514) reflecting boundaries for parameter proposals. Now
    % all proposals will be run.
    for cc=1:MCMCPar.n
        while any(X_prop(:,cc)<MinMax(:,1))>0 | any(X_prop(:,cc)>MinMax(:,2))>0 %if any are outside
            too_low = find(X_prop(:,cc)<MinMax(:,1));
            too_high = find(X_prop(:,cc)>MinMax(:,2));
            diff_low = -X_prop(too_low,cc)+MinMax(too_low,1);
            diff_high= -X_prop(too_high,cc)+MinMax(too_high,2);
            X_prop(too_low,cc)  = MinMax(too_low,1)+diff_low;
            X_prop(too_high,cc) = MinMax(too_high,2)+diff_high;
        end
    end
    
    doruns = 1:MCMCPar.n;
%     doruns(ignorenow)=[]; %Removed since boundaries are reflecting.
%     JOS060514
 
    % parellel bits
   
    X_actual = 10.^(X_prop(:,doruns)); 
    nlogl_tmp = NaN*size(nlogl_new); 
    SS_tmp = NaN*size(SS_new); 
    
   % for iter_parfor = 1:num_cores:numel(doruns)
   fprintf('Iteration %d, %d chains \n', tt, numel(doruns)) 
   
    for dd=1:numel(doruns) 
        %parfor dd = iter_parfor:(iter_parfor+num_cores-1) 
        [nlogl_tmp(dd,:),SS_tmp(dd,:),Nobs] = MyLake_mcmc_run(X_actual(:,dd),partix,fixix,fixval,K_values_sediments,K_values_lake) 
        %pause(2)
        %disp('cooling down')
    end
    
    %total_runs = total_runs * num_cores;
  
    
    %end
    nlogl_new(doruns,:) = nlogl_tmp(1:numel(doruns),:);
    SS_new(doruns,:)    = SS_tmp(1:numel(doruns),:); 

    
    alphas(doruns) = min(1,exp(-sum(nlogl_new(doruns,:),2)+sum(nlogl(doruns,:),2)));
    swap = find(alphas>rand(MCMCPar.n,1));
    Accept(tt) = numel(swap)/MCMCPar.n;
    X_now(:,swap) = X_prop(:,swap);
    nlogl(swap,:) = nlogl_new(swap,:);
    SS_Hist(:,:,tt) = SS_new;
    nLogHist(:,tt) = sum(nlogl,2);
    Chains(:,:,tt) = X_now;
    if tt/MCMCPar.checkGR==floor(tt/MCMCPar.checkGR)
        % if checking Gelman Rubin convergence
        [R_stat] = GelmanRubin_v2(Chains,MCMCPar,tt);
        GR_Hist(:,gr_tix) = R_stat';
        gr_tix = 1;
    end
	
	%Essentially the algorithm fundoutliers_RMV_v1.m compare the history of the likelihoods of the chains, 
	%and removes the very unlikely ones and moves those chains to more likely parts of 
	%the parameter space. 
	
	if tt<MCMCPar.burnin & (floor(tt/MCMCPar.checkoutlier)==tt/MCMCPar.checkoutlier) 
		[nLogHist,Chains,SS_Hist,SS,nlogl,outlierc] = findoutliers_RMV_v1(nLogHist,Chains,SS_Hist,SS,tt,nlogl,outlierc); 
	end 
	    
    if tt<10
        figure(1)
        tmp = squeeze(Chains(1,:,1:tt));
        subplot(2,2,1),plot(tmp')
        tmp = squeeze(Chains(2,:,1:tt));
        subplot(2,2,2),plot(tmp')
        tmp = squeeze(Chains(3,:,1:tt));
        subplot(2,2,3),plot(tmp')
        subplot(2,2,4),plot(Accept(1:tt))
    elseif tt/plt ==floor(tt/plt)
        disp('backing up workspace ...')
        save par_dream_temp.mat  
        figure(1)
        tmp = squeeze(Chains(1,:,1:tt));
        subplot(2,2,1),plot(tmp')
        tmp = squeeze(Chains(2,:,1:tt));
        subplot(2,2,2),plot(tmp')
        tmp = squeeze(Chains(3,:,1:tt));
        subplot(2,2,3),plot(tmp')
        subplot(2,2,4),plot(Accept(1:tt))
       
    end
end

toc
%%
% So we need to
% 1. detail the first run and collect nlogl SS, and sigmas
% 2. generate the main loop. inside this loop;
%   - generate proposal values.
%   - run model
%   - do MH
%   - storing
%   - testing GR checks
%   - outlier checks

