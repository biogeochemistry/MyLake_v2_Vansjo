
function [par_file,dat_file,lfert_file,obs_file,ssf_file,out_file]= fn_INCA_generate_filenames(runID)
%% will generate filename set with suffixes 'reachID' and 'runID'. They can be omitted in the function call. 

par_file =      (['vansjo' runID '.par']);
dat_file =      (['vansjo' runID '.dat']);
lfert_file =    (['vansjo' runID '.lfd']);
obs_file =      (['vansjo' runID '.obs']);
ssf_file =      (['vansjo' runID '.ssf']);
out_file =      (['vansjo' runID '.dsd']);

end