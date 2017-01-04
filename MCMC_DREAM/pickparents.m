function [r1s r2s] = pickparents(no_chains,delta)
% [r1s r2s] = pickparents(no_chains,delta)
%This function picks delta numbers of other chains to use in the DREAM
%algorithm in producing candidate points. There are no_chains and for each
%chain this function draws delta number of pairs.

for cc=1:no_chains;
    notthis = 1:no_chains;
    notthis(cc) = [];
    for nn=1:delta
        % QUESTION: should it never be the same TWO or the same across
        % all? Now across all
        [a] = randi(length(notthis),1); %r1
        r1s(cc,nn) = notthis(a);
        notthis(a) = []; %removing this
        [b] = randi(length(notthis),1); %r2
        r2s(cc,nn) = notthis(b);
        notthis(b) = [];
    end
end