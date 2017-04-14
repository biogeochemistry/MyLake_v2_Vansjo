function [MyLake_results, Sediment_results] = fn_MyL_application(m_start,m_stop, K_sediments, K_lake, use_INCA, run_INCA, run_ID, clim_ID)
global sed_par_file lake_par_file Eevapor
% This is the main MyLake application configuration file. INCA is a switch
% It is made to run a after the parameter are set by Set_Prior

Eevapor=0;    disp('init ...');

calibration_k_values = [(1:52)',cell2mat(K_sediments(:,1)) ]; % writing sediments parameters file

%% generates unique files

sed_par_file = tempname;
lake_par_file = tempname;

dlmwrite(sed_par_file, calibration_k_values,'delimiter','\t');

%% writing lake parameter file
f = fopen('IO/vansjo_para.txt');
garbage = fgetl(f); % file get line
garbage = fgetl(f); % file get line
data_lake = textscan(f, '%s%f%f%f%s', 64, 'Delimiter', '\t');
fclose(f); % the parameter line (xx,1) + 2 lines gives the location of the paramter in the input txt file.
% array position + 2 = input file line

data_lake{1, 2}(1,1) = K_lake{1}; % I_scDOC
data_lake{1, 2}(2,1) = K_lake{2}; % qy_DOC
data_lake{1, 2}(3,1) = K_lake{3}; % oc_DOC
data_lake{1, 2}(4,1) = K_lake{4}; % k_bod
data_lake{1, 2}(5,1) = K_lake{5}; % theta_bod
data_lake{1, 2}(6,1) = K_lake{6}; % theta_bod_ice
data_lake{1, 2}(7,1) = K_lake{7};
data_lake{1, 2}(8,1) = K_lake{8};
data_lake{1, 2}(9,1) = K_lake{9};
data_lake{1, 2}(10,1) = K_lake{10}; % BOD7  on line 55 of the txt file
data_lake{1, 2}(11,1) = K_lake{11}; % I_scO
data_lake{1, 2}(12,1) = K_lake{12}; % Q10
data_lake{1, 2}(13,1) = K_lake{13}; % wc_factor
data_lake{1, 2}(14,1) = K_lake{14}; % 
data_lake{1, 2}(15,1) = K_lake{15}; % 
data_lake{1, 2}(16,1) = K_lake{16}; % 
data_lake{1, 2}(17,1) = K_lake{17}; % 
data_lake{1, 2}(18,1) = K_lake{18}; % 
data_lake{1, 2}(19,1) = K_lake{19}; % 
data_lake{1, 2}(20,1) = K_lake{20}; % 
data_lake{1, 2}(21,1) = K_lake{21}; % 
data_lake{1, 2}(22,1) = K_lake{22}; % 
data_lake{1, 2}(23,1) = K_lake{23}; % 
data_lake{1, 2}(24,1) = K_lake{24}; % 
data_lake{1, 2}(25,1) = K_lake{25}; % 
data_lake{1, 2}(26,1) = K_lake{26}; % 
data_lake{1, 2}(27,1) = K_lake{27}; % 
data_lake{1, 2}(28,1) = K_lake{28}; % 
data_lake{1, 2}(29,1) = K_lake{29}; % 
data_lake{1, 2}(30,1) = K_lake{30}; % 
data_lake{1, 2}(31,1) = K_lake{31}; % 
data_lake{1, 2}(32,1) = K_lake{32}; % 
data_lake{1, 2}(33,1) = K_lake{33}; % 
data_lake{1, 2}(34,1) = K_lake{34}; % 
data_lake{1, 2}(35,1) = K_lake{35}; % 
data_lake{1, 2}(36,1) = K_lake{36}; % 
data_lake{1, 2}(37,1) = K_lake{37}; % 
data_lake{1, 2}(38,1) = K_lake{38}; % 
data_lake{1, 2}(39,1) = K_lake{39}; % 
data_lake{1, 2}(40,1) = K_lake{40}; % 
data_lake{1, 2}(41,1) = K_lake{41}; % 
data_lake{1, 2}(42,1) = K_lake{42}; % 
data_lake{1, 2}(43,1) = K_lake{43}; % 
data_lake{1, 2}(44,1) = K_lake{44}; % 
data_lake{1, 2}(45,1) = K_lake{45}; % 
data_lake{1, 2}(46,1) = K_lake{46}; % 
data_lake{1, 2}(47,1) = K_lake{47}; % 
data_lake{1, 2}(48,1) = K_lake{48}; % 
data_lake{1, 2}(49,1) = K_lake{49}; % 
data_lake{1, 2}(50,1) = K_lake{50}; % 
data_lake{1, 2}(51,1) = K_lake{51}; % 
data_lake{1, 2}(52,1) = K_lake{52}; % 
data_lake{1, 2}(53,1) = K_lake{53}; % 
data_lake{1, 2}(54,1) = K_lake{54}; % 
data_lake{1, 2}(55,1) = K_lake{55}; % 
data_lake{1, 2}(56,1) = K_lake{56}; % 
data_lake{1, 2}(57,1) = K_lake{57}; % 
data_lake{1, 2}(58,1) = K_lake{58}; % 
data_lake{1, 2}(59,1) = K_lake{59}; % 
data_lake{1, 2}(60,1) = K_lake{60}; % 
data_lake{1, 2}(61,1) = K_lake{61}; % 
data_lake{1, 2}(62,1) = K_lake{62}; % 
data_lake{1, 2}(63,1) = K_lake{63}; % 
data_lake{1, 2}(64,1) = K_lake{64}; % 


fid=fopen(lake_par_file,'wt');
fprintf(fid,'\n\n');
dlmwrite(lake_par_file, [[1:64]',data_lake{2},data_lake{3},data_lake{4},(1:64)'],'delimiter','\t','-append'); % 1:64 is the length of the parameter file.
fclose(fid);

%% Specific MyLake application

% warning('off', 'all')
lake='Vansjo';
year=1983;
dt = 1.0;

%% If you want MyLake to read INCA outputs use_INCA
% (1) Reading exisitng INCA outputs to prepare MyLake input
% (0) You are running only MyLake, inputs alrady exist in the IO folder
if use_INCA == 1;
   [store_INCAP_input, vanem_INCAP_input, INCA_QC]=fn_INCA_MyL(run_INCA, run_ID, clim_ID, m_start, m_stop);
else
    INCA_QC = 0;
end


% store_INCAP_input
% vanem_INCAP_input


%# ############ This is Vansjø Storefj ##############

parafile=lake_par_file;
initfile='IO/mylake_initial_concentrations.txt';


if use_INCA == 0
    inputfile='IO/store_INCAP_input_baseline_mod.txt';
    disp('Using existing input')
elseif use_INCA == 1
    inputfile = store_INCAP_input; % setting use_INCA to 2 will look for store_INCAP_input
    disp('Using INCA output')
elseif ischar(use_INCA);
    inputfile=use_INCA;
    disp('Using response surfaces array')
end


Deposition = 0;

disp('Storefjorden ...')

[MyLake_results_basin1, sediment_data_basin1] ...
    = solvemodel_v2(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake');

%% Old KOJI stuff
% date_temp = datevec(tt);
% date_temp = date_temp(:,1:3);
% %csvwrite('temp/dateinformation.csv', [m_start; m_stop]);
% %% Write lake physics ... disable to increase speed of MCMC run
% % dlmwrite('basin1Tzt.csv', [date_temp Tzt'],'delimiter', ',', 'precision', '%5.2f');
% % dlmwrite('basin1Qst.csv', [date_temp Qst'],'delimiter', ',', 'precision', '%5.2f');
% % dlmwrite('basin1Wt.csv', [date_temp Wt], 'delimiter', ',', 'precision', '%5.2f');
% % dlmwrite('basin1His.csv', [date_temp His'],'delimiter', ',', 'precision', '%5.4f');
% % dlmwrite('basin1freezingdates.csv', DoF','delimiter', ',');
% % dlmwrite('basin1meltingdates.csv', DoM','delimiter', ',');
%
% % Qst : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
% #################################################################

% surfacearea = Az(1); % m2
% precipvolume = surfacearea * Wt(:, 7) / 1000; % m3 day-1
% runoffintolake = Phys_par(16) * Inflw(:, 1); % I_scV should be 1'is input scaling 1?'
% outflow = precipvolume + runoffintolake; %% dd x 1
% outflowTemp = Tzt(1, :)';
% outflowC = Czt(1, :)';
% outflowS = Szt(1, :)';
% outflowTP = Czt(1, :)' + Pzt(1, :)' + Chlzt(1, :)' + PPzt(1, :)' + DOPzt(1, :)' ;
% outflowDOP = DOPzt(1, :)';
% outflowChl = Chlzt(1, :)';
% outflowDOC = DOCzt(1, :)';
% outflowDIC = DOCzt(1, :)'; %dummy for MyLake TSA
% outflowO = DOCzt(1, :)'; %dummy for MyLake TSA
% outflowDIC = DICzt(1, :)';
% outflowO = O2zt(1, :)';
% outflowNO3zt = NO3zt(1,:)';
% outflowNH4zt = NH4zt(1,:)';
% outflowSO4zt = SO4zt(1,:)';
% outflowHSzt = HSzt(1,:)';
% outflowH2Szt = H2Szt(1,:)';
% outflowFe2zt = Fe2zt(1,:)';
% outflowCa2zt = Ca2zt(1,:)';
% outflowpHzt = pHzt(1,:)';
% outflowCH4zt = CH4zt(1,:)';
% outflowFe3zt = Fe3zt(1,:)';
% outflowAl3zt = Al3zt(1,:)';
% outflowSiO4zt = SiO4zt(1,:)';
% outflowSiO2zt = SiO2zt(1,:)';
% outflowdiatomzt = diatomzt(1,:)';


% %# ############ This is Vansjø Vanemfj. ##############
% if isnumeric(use_INCA) % to avoid running two basins in case of RS analysis. 
    
%     if use_INCA == 0
%         land_to_vanem = 'IO/vanem_INCAP_input_baseline_mod.txt'
%     else
%         land_to_vanem = vanem_INCAP_input;  % created above by calling fn_INCA_MyL.m
%     end
    
    
%     store_to_vanem = [outflow outflowTemp outflowC outflowS outflowTP outflowDOP outflowChl outflowDOC outflowDIC outflowO outflowDIC outflowO outflowNO3zt outflowNH4zt outflowSO4zt outflowHSzt outflowH2Szt outflowFe2zt outflowCa2zt outflowpHzt outflowCH4zt outflowFe3zt outflowAl3zt outflowSiO4zt outflowSiO2zt outflowdiatomzt];
   
%     Q_lake = outflow;
    
%     vanem_input = tempname;
%     merge_l_b_inputs(land_to_vanem,store_to_vanem,vanem_input)
    
%     %parafile='k_values_lake.txt';
%     parafile = lake_par_file;
%     initfile='IO/vanem_init.txt';
%     inputfile = vanem_input
    
%     % note: I removed the DIC/O2 bits here ... take them again from Langtjern
%     % app when migrating to Mylake DOCOMO
    
    
    
%     disp('Vanemfjorden ...')
    
%     [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DICz,In_DOCz,In_TPz_sed,In_Chlz_sed,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_FIM,Ice0,Wt,Inflw,...
%     Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
%         = modelinputs_v2(m_start,m_stop, initfile, 'lake', inputfile, 'timeseries', parafile, 'lake', dt);
    
%     [zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,DICzt,CO2zt,O2zt,NO3zt,NH4zt,SO4zt,HSzt,H2Szt,Fe2zt,Ca2zt,pHzt,CH4zt,Fe3zt,Al3zt,SiO4zt,SiO2zt,diatomzt,O2_sat_relt,O2_sat_abst,BODzt,Qzt_sed,lambdazt,...
%     P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,surfaceflux,O2fluxt,CO2_eqt,K0t,O2_eqt,K0_O2t,...
%     CO2_ppmt,dO2Chlt,dO2BODt,testi1t,testi2t,testi3t,...
%     MyLake_results_basin2, sediment_data_basin2] ...
%         = solvemodel_v2(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake',In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_TPz_sed,In_Chlz_sed,In_FIM, ...
%     Ice0,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names, Bio_par,Bio_par_range,Bio_par_names, Deposition);
   
%     % delete (vanem_input)
%     disp('Cleanup ... done.')
% end
 

% MyLake_results = [MyLake_results_basin1, MyLake_results_basin2];
MyLake_results = [MyLake_results_basin1];
% Sediment_results = [sediment_data_basin1, sediment_data_basin2];
Sediment_results = [sediment_data_basin1];


%% cleaning
fclose('all');
% delete (sed_par_file)
% delete (lake_par_file)

end
