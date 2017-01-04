function [L,X_prop,WhichM] = generateProposalv4(MCMCPar,X,L,tt,pCR)
% .m file to generate proposal values according to the MCMC-dream
% algorithm, used for the FW estimation in AMAP.
% v4 - fixed the gamma=1 problematic and allows now for MCMCPar.delta = 0.
% - when gamma==1 only compare with one pair and replace all,
%
% MCMCPar structure uses entries
% .delta
% .n - no chains
% .nparam - no parameters (excluding gibbs sampled ones)
% .b_star - randn(0,b_start) for ergodicity
% .b - uniform scaling the jumps
% .CRupdate - Yes or No
% .nCR
% .gammaisone
% CR values need to be defined
% pCR list of crossover probs e.g. [1/3 2/3 3/3]
% X have dimensions (nopar,nochain)
% L are for updating the CR values.
% tt is iteration number for gammaisone eval.


if MCMCPar.delta >0;
    [r1s,r2s] = pickchains(size(X,2),MCMCPar.delta);
end
epss = randn(MCMCPar.nparam,MCMCPar.n).*(MCMCPar.b_star^2);
%random('Normal',0,MCMCPar.b_star,MCMCPar.nparam,MCMCPar.n); %moving the draw outside the loop for efficiency.
es = -MCMCPar.b + 2*MCMCPar.b*rand(MCMCPar.nparam,MCMCPar.n); %e from Vrugt
% Reduction in time by 12/89 s --> 10 %.
for cc=1:MCMCPar.n % for each chain.
    if strcmp(MCMCPar.CRupdate,'Yes') %if updating the crossover probs
        mrn = mnrnd(1,pCR);
        L = L +mrn;
        WhichM(cc) = find(mrn);
        CR_now = find(mrn)/MCMCPar.nCR; %which crossover value now.
    else
        CR_now = 1/3; %default.
        L = L;
        WhichM = 1; %for stability.
    end
    
    %     x_p = NaN*ones(MCMCPar.nparam,1); %put proposals here
    
    if (floor(tt/MCMCPar.gammaisone))==tt/MCMCPar.gammaisone %if every xth generation
        % This is randomized in Vrugts code, with rand<P(replace)
        gamma = 1;
        if MCMCPar.delta>0 % if we do look at other chains
            
            delta_now = gamma*(sum(X(:,r1s(cc,1)),2)-sum(X(:,r2s(cc,1)),2)); %only sampling one comparison
        else
            delta_now = 0;
        end
        X_prop(:,cc) = (X(:,cc) + (1 + es(:,cc)).*delta_now+epss(:,cc));
        % Now swapping all.
        
        %UPDATED BY JOS 270913, seeing how Vrugt does in his code.
    else
        replace = binornd(1,1-CR_now,MCMCPar.nparam,1);
        if sum(replace) ==length(replace) %if all are replaced, at least one sohuld be changed
            replace(randi(length(replace)))=0;
        end
        d_dash = MCMCPar.nparam-sum(replace);
        
        gamma = 2.38/(sqrt(2*MCMCPar.delta*d_dash));
        if MCMCPar.delta>0
            usedelta = randi(MCMCPar.delta);
            delta_now = gamma*(sum(X(:,r1s(cc,1:usedelta)),2)-sum(X(:,r2s(cc,1:usedelta)),2));
        else
            delta_now = 0;
        end
        X_prop(:,cc) = (X(:,cc) + (1 + es(:,cc)).*delta_now+epss(:,cc)).*(1-replace) + ...
            replace.*X(:,cc);
        % putting some back here, according to CR value and draw.
        
    end
    
    % ACtually in Vrugts code, when gamma is set to 1 he only looks at the
    % difference between two chains. This makes sense, since when gamma~=0
    % it is the sum of the difference of many chains.
    
    
    
    
    % Taken out of loop
    %     for pp=1:MCMCPar.nparam %for each parameter
    %         e = -MCMCPar.b + 2*MCMCPar.b*rand; %e from Vrugt
    %         eps = random('Normal',0,MCMCPar.b_star); %now not scaled by min-max as in INCA application
    
    %         x_p(pp) = X(pp,cc) + (1+es(cc,pp))*gamma*(sum(X(pp,r1s(cc,1:usedelta)))-sum(X(pp,r2s(cc,1:usedelta)))) + epss(cc,pp);
    %         X_prop(pp,cc) = (replace(pp)).*X(pp,cc) + (1-replace(pp))*x_p(pp);
    
    %     end
    % Can we do this outside of the loop?
    % size of arrays
    % replace ( 1: pp)
    % es ( 1: cc, 1:pp)
    % epss (1:cc, 1:pp)
    % X (1:pp,1:cc))
    %
end