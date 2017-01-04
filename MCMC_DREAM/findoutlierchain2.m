function [Param_Track,X,outlier,Hist_logp,SS_conc_old] = findoutlierchain2(X,Hist_logp,outlier,Param_Track,SS_conc_old,MCMCPar,tt)
%[Param_Track,X,outlier,Hist_logp] = findoutlierchainv2(X,Hist_logp,outlier,Param_Track,MCMCPar,tt)
%Finding outlier chains using the Inter-Quartile-Range statistics.
%Param_Track is now a cell array Param_Tracl{chain}(parameter,time)

%was (parameter, chain, time)
%Hist_logp   is (chain, time)
idend = tt; %length of simulation currently
idst  = floor(idend/2);    %start of chain for which the stats are calculated

mean_hist_logp = nanmean(Hist_logp(:,idst:idend),2); %Mean log density

% Derive the upper and lower quantile of the data
Q1 = prctile(mean_hist_logp,75); 
Q3 = prctile(mean_hist_logp,25);
% Derive the Inter quartile range
IQR = Q1 - Q3;
% Compute the upper range -- to detect outliers
UpperRange = Q3 - 2 * IQR;
% See whether there are any outlier chains
chain_id = find(mean_hist_logp < UpperRange);
Nid = size(chain_id,1);

if (Nid > 0),
    % Loop over each outlier chain
    for qq = 1:Nid,
        % Finding the best chain.
        r_idx = find(mean_hist_logp==max(mean_hist_logp)); r_idx = r_idx(1);
        % Added -- update hist_logp -- chain will not be considered as an outlier chain then
        Hist_logp(chain_id(qq),1:end) = Hist_logp(r_idx,1:end);
        % Jump outlier chain to r_idx -- Sequences
        %Param_Track(1:MCMCPar.nparam,chain_id(qq),1:idend) = Param_Track(1:MCMCPar.nparam,r_idx,1:idend); %substituting the whole param tracking for future outlier calcs.
        %The call above uses too much memory, will it be easier to first
        %copy the r_idx params and then put them into the Param_Track,
        %since matlab then need not make a full copy of the whole
        %Param_Track?
        Param_Track{chain_id(qq)} = Param_Track{r_idx}; %swapping chains.
        SS_conc_old(chain_id(qq)) = SS_conc_old(chain_id(qq)); %swapping stored sums of squares.
        %Param_Track(1:MCMCPar.nparam,chain_id(qq),1:idend) = tmp;
        
        
        % Old
        %         tmp = Param_Track(1:MCMCPar.nparam,r_idx,1:idend);
%         Param_Track(1:MCMCPar.nparam,chain_id(qq),1:idend) = tmp;
%         
        
        %X(r_idx,1:MCMCPar.nparam); %Putting in the Position at the current timestep. NOTE THIS IS DIFFERENT FROM VRUGTS CODE
        % Jump outlier chain to r_idx -- X
        X(1:MCMCPar.nparam,chain_id(qq)) = X(1:MCMCPar.nparam,r_idx);
        % Add to chainoutlier
        outlier = [outlier ; tt chain_id(qq)];
    end;
else
    % Do nothing
end