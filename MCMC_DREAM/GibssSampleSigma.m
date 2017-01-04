function [sig] = GibssSampleSigma(n0,ss0,ss,Nobs)
% This uses n0 (default 1) and ss0 (given initial error variance) and
% actual sumsofsquars (ss) together with number of datapoints used to get
% the ss. This is one dimensional, but accepts vectors/columns.

% Cut from INCA application, Gibbs sampling of error terms.
% n0 = [10 10 10];
% SS0 = [50 150 50];
sig = NaN*ones(size(n0));

for ii= 1:numel(n0)
    %     c2tmp = chi2rnd(Nobs(ii));
    %Sigmas(cc,ii) = sqrt(1/c2tmp*SS_conc_old(cc,ii)); %Following Vrugt et al 'Equifinality ...'
    %NOTE now using sqrt since Sigmas is sigma and not sigma^2.
    
    %FROM TSA These are identical.
    %         Sigmas(2)=1./sqrt(gamrnd((0.001+Nobs_conc(2)/2), 1./(0.001+0.5*SS_conc_new(2))));
    %                 Sigmas(1)=1./sqrt(gamrnd((0.001+Nobs_conc(1)/2), 1./(0.001+0.5*SS_conc_new(1))));
    %TESTING TESTING TESTING 031011
    sig(ii)=1./sqrt(gamrnd(((n0(ii)+Nobs(ii))/2), 1./((ss0(ii)+ss(ii))/2)));
end
