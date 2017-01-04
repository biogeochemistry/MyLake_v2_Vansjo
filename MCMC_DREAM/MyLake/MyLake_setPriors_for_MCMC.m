% Generating Prior structure. For this application we will have an entry
% into the Prior even for parameters not varied. Put in a switch for that.

MaxMinSed = ...
    [1 1; ... 
    0.01 0.01; ...
    0.0123 0.0123;...
    0.01 0.01;...
     3.92  3.92;...
    2415 2415; ...
    0.0293 0.0293;...
    0.001 0.001; ...
    0.1 0.1;...
    1e-2 1e-2;...
    0.1 0.1;...
    0.1 0.1; ...
    0.1 0.1;...
    100 100; ...
    1 1; ...
    0.1 0.1;...
    0.1 0.1; ...
    3.17 3.17;...
    0.1 0.1;...
    1.35 1.35;...
    1.35 1.35; ...
    6500 6500;...
    0.1 0.1; ...
    2500 2500; ...
    0.001 0.001;...
    21.3 21.3; ...
    0.37 0.37;...
    3e-6 3e-6; ...
    8e-2 8e-2; ...
    1000 1000;...
    0.001 0.001;...
    0 0];


K_values_sediments ...
    = {0,   'K_OM1'; 	% optimized
    0,      'K_OM2'; 		% optimized
    0.0123, 'Km_O2'; 	% optimized 
    0.01,   'Km_NO3'; 	% not varied
    3.92,   'Km_Fe(OH)3';% optimized
    2415,   'Km_FeOOH'; 	% optimized
    0.0293, 'Km_SO4'; 	% not varied
    0.001,  'Km_O2ao'; 	% optimized
    0.1,    'Km_NH4ao'; 	% not varied
    0.0219, 'Kin_O2'; 	% 1e-2, 1
    0.1,    'Kin_NO3'; 	% not varied
    0.1,    'Kin_FeOH3'; 	% not varied
    0.1,    'Kin_FeOOH'; 	% not varied
    100,    'K_NH4ox'; 	% not varied
    8.7e4,  'K_Feox'; 	% optimized
    0.1,    'K_Sdis'; 	% not varied
    0.1,    'K_Spre'; 	% not varied
    3.17,   'K_FeS2pre'; 	% not varied
    0.1,    'K_AlOH3'; 	% not varied
    1.35,   'K_P_sorb_a'; 	% not varied
    1.35,   'K_P_sorb_b'; 	% not varied
    6500,   'K_rhom'; 	% not varied
    0.1,    'K_tS_Fe';  	% not varied
    2500,   'K_Fe_S'; 	% not varied
    0.001,  'K_Fe_dis'; 	% not varied
    21.3,   'K_Fe_pre'; 	% not varied
    0.37,   'K_apa'; 		% not varied
    3e-6,   'Kapa'; 	% not varied
    8e-2,   'K_sorg'; 	% not varied
    1000,   'K_tsox'; 	% not varid
    0.001,  'K_FeS_FeS2'; 	% not varied
    0,      'accel'};		% , 0.01-50 #1

MaxMinLake =     ... 
    [0 2; ...       % I_scT
    1 1; ...        % I_scDOC
    0.2 0.2; ...    % w_s
    0.2 0.2; ...    % w_chl
    0.7 1.3; ...        % Y_cp
    0.1 0.4; ...    % m_twty
    0.75 3.0; ...    % g_twty
    1e-4 4e-4; ...  % k_sed_twty
    0 0; ...        % k_dop_twty
    0.1 0.4; ...    % P_half
    0.01 0.01; ...  % oc_DOC
    0.1 0.1; ...    % qy_DOC
    0.1 0.1; ...    % k_BOD
    1.05 1.05; ...  % theta_BOD
    1.13 1.13; ...  % theta_BOD_ice
    4 4; ...        % theta_T
    1 1];           % I_scO

K_values_lake =      ...
   {0,      'I_scT'; 	% 
    1,      'I_scDOC'; 	% 
    0.2,    'w_s'; 		% 
    0.2,    'w_chl';    % 
    1,      'Y_cp'      %
    0.2,    'm_twty'    %
    1.5     'g_twty'     %
    2e-4    'k_sed_twty' %
    0       'k_dop_twty' %    
    0.2     'P_half'     %
    0.01,   'oc_DOC'; 	 % 
    0.1,    'qy_DOC'; 	 % 
    0.1,    'k_BOD'; 	 % 
    1.05,   'theta_bod';	 % 
    1.13,   'theta_bod_ice'; 	% 
    4,      'theta_T';     % 
    1,      'I_scO'}; 	% 

for ii=1:32 %for sediment priors
    Priors(ii).name = K_values_sediments(ii,2);
    Priors(ii).min  = MaxMinSed(ii,1);
    Priors(ii).max  = MaxMinSed(ii,2);
    Priors(ii).default = K_values_sediments(ii,1);
end

for ii=1:17 %for lake parameters
    Priors(32+ii).name = K_values_lake(ii,2);
    Priors(32+ii).min = MaxMinLake(ii,1);
    Priors(32+ii).max = MaxMinLake(ii,2);
    Priors(32+ii).default = K_values_lake(ii,1);
end

