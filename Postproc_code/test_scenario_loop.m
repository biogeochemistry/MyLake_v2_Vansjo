state_var = {...
    'dummy',...
    'TP',...
    'Chl',...
    'PO4',...
    'Ppart',...
    'T'};

% loop to nowhere , eventually loop thorugh all vars
for cur_sv = 2:length(state_var);
    for current_run = 1:no_runs

Scenario = {...
'NO_Vansjo_Hist_',...
'NO_Vansjo_Base_'...
'NO_Vansjo_Hist_M2_'...
'NO_Vansjo_Hist_M3_'...
'NO_Vansjo_G4_',...
'NO_Vansjo_I4_',...
'NO_Vansjo_G8_',...
'NO_Vansjo_I8_',...
'NO_Vansjo_G8Tech_',...
'NO_Vansjo_G8Cons_',...
'NO_Vansjo_G8Frag_',...
'NO_Vansjo_I8Tech_',...
'NO_Vansjo_I4Cons_',...
'NO_Vansjo_I8Frag_',};
        
run_ID = strcat(Scenario{current_run},state_var(cur_sv)) %  CALIBRATION RUN

    end
end
