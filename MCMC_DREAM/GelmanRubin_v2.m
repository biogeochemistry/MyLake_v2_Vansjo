function [R_stat] = GelmanRubin_v2(Chains,MCMCPar,tt)
% Chains is the array of parameters with dimensions (pars, chain, iteration)
%
idend = tt;
idst  = ceil(tt/2); %only using the last half of the series.

for cc=1:MCMCPar.n %for each chain
    for pp=1:MCMCPar.nparam %for each parameter
    VarWithin(cc,pp) = var(Chains(pp,cc,idst:idend)); %
%     var(Param_Track{cc}(pp,idst:idend)); %
    MeanWithin(cc,pp) = mean(Chains(pp,cc,idst:idend));
    %mean(Param_Track{cc}(pp,idst:idend));
    end
end

for pp=1:MCMCPar.nparam %For each parameter
    VarBetween(pp) = var(MeanWithin(:,pp)); %variance in means
    B(pp) = MCMCPar.n*VarBetween(pp);
end

W = mean(VarWithin); %mean variance across chains.

sigma2 = ((MCMCPar.n-1)/(MCMCPar.n))*W + (1/MCMCPar.n)*B;
R_stat = sqrt(((MCMCPar.nparam+1)/MCMCPar.nparam)*sigma2./W - (MCMCPar.n-1)/MCMCPar.nparam/MCMCPar.n);
