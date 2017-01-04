function [nLogHist,Chains,SS_Hist,SS,nlogl,outlierc] = findoutliers_RMV_v1(nLogHist,Chains,SS_Hist,SS,tt,nlogl,outlierc)
% A quick little outlier detection for RMC's MyLake. This will compute the
% average logarithm of the likelihoods the last half of the simulation, and
% move chains that are 2 times the interquartile range from the least
% likely to the chain with the most likely position. 

% 1. check for outliers
% 2. Identify/index the outlier AND the most likely.
% 3. Swap all histories and current pars with most likely.

% I: [For 1&2] nlogl(chain, type_obs)
% [For 3] nLogHist(chain,iter), Chains(par,chain,iter),
% SS_Hist (chain, type_obs, iter), SS(chain,typeobs),tt
% [For control]: outlier_counter
% O: all of the above + outlier-counter.

% 1.
%Finding outlier chains using the Inter-Quartile-Range statistics.
%Hist_logp   is (chain, time)
idend = tt;                 % length of simulation currently
idst  = floor(idend/2);     % start of chain for which the stats are calculated
mean_hist_logp = nanmean(nLogHist(:,idst:idend),2);
% Mean negative log likelihood.

% Derive the upper and lower quantile of the data
Q_u = prctile(mean_hist_logp,75);
Q_l = prctile(mean_hist_logp,25);
% Derive the Inter quartile range
IQR = Q_u - Q_l;
% Compute the upper range -- to detect outliers
UpperRange = Q_u + 2 * IQR;

Out = find(mean_hist_logp>UpperRange);
% These are outliers and will be moved to best chain.
if numel(Out)>0 %if any outliers
    [~,best] = min(mean_hist_logp);
    for ii=1:numel(Out) %for each of the outliers detected
        % Changing the histories
        nLogHist(Out(ii),:)     = nLogHist(best,:);
        Chains(:,Out(ii),:)     = Chains(:,best,:);
        SS_Hist(Out(ii),:,:)    = SS_Hist(best,:,:);
        SS(Out(ii),:)           = SS(best,:);
        nlogl(Out(ii),:)        = nlogl(best,:);
        outlierc                = outlierc+1;
    end
    
end

%
