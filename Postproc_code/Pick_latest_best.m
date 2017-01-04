[~,best]=min(sum(nlogl,2));
Temp = Chains(:,best,end);
Param_now = 10 .^ Temp