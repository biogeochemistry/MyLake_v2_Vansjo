tic
% parpool
gaoptions = optimoptions('ga','UseParallel',true);
x0 = [2500, 8000, 0.02, 0.2, 1.5, 0.2, 0.02  0.2, 1.5, 0.2, 2]
x = fmincon(@opt_fun,x0)


% rng default % to get the same evaluations as the previous run
% % if gaAvailable
%     % gaoptions = optimoptions(gaoptions,'UseParallel',true);
%     startTime = tic;
%     gasol = ga(@opt_fun,11,[],[],[],[],[],[],[],gaoptions);
%     time_ga_parallel = toc(startTime);
%     fprintf('Parallel GA optimization takes %g seconds.\n',time_ga_parallel);
% % end
