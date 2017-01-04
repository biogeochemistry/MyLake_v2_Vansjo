function [outlier] = findoutlierchain_v3(Hist_logp,tt)
% This is a new function which does not change anything in the arrays (as
% was done pre180214, but only gives an output (outlier) which is:
% [tt outlier_inx best_inx]
% This means that all chains and SS and whatnot that needs to be changed
% when removing an outlier can be done using the values here.
% This version (v3) also moves chains that are at the minimum log
% likelihood, defined as <log(realmin*2), the smallest number for matlab.

% Hist_logp is an array with stored log likelihood values (iteration by chain)
% Note that this function ignores NaNs and will therefore not move any
% chains with NaN as loglikelihood
outlier = [];

%Finding outlier chains using the Inter-Quartile-Range statistics.
%Hist_logp   is (chain, time)
idend = tt;                 % length of simulation currently
idst  = floor(idend/2);     % start of chain for which the stats are calculated
mean_hist_logp = nanmean(Hist_logp(idst:idend,:),1); %Mean log density

% Derive the upper and lower quantile of the data
Q1 = prctile(mean_hist_logp,75); 
Q3 = prctile(mean_hist_logp,25);
% Derive the Inter quartile range
IQR = Q1 - Q3;
% Compute the upper range -- to detect outliers
UpperRange = Q3 - 2 * IQR;
% See whether there are any outlier chains
chain_id = find(mean_hist_logp < UpperRange);
% For version 3
chain_id = [find(mean_hist_logp < UpperRange) find(mean_hist_logp<log(realmin*2))]
Nid = size(chain_id,2);

if (Nid > 0),
    % Loop over each outlier chain
    for cc = 1:Nid,
        % Finding the best chain.
        r_idx = find(mean_hist_logp==max(mean_hist_logp)); 
        r_idx = r_idx(1);
        % Added -- update hist_logp -- chain will not be considered as an outlier chain then
%         Hist_logp(chain_id(cc),1:end) = Hist_logp(r_idx,1:end);
        % Jump outlier chain to r_idx -- Sequences
        %Param_Track(1:MCMCPar.nparam,chain_id(qq),1:idend) = Param_Track(1:MCMCPar.nparam,r_idx,1:idend); %substituting the whole param tracking for future outlier calcs.
        %The call above uses too much memory, will it be easier to first
        %copy the r_idx params and then put them into the Param_Track,
        %since matlab then need not make a full copy of the whole
        outlier = [outlier ; tt chain_id(cc) r_idx];
    end;
else
    % Do nothing
end