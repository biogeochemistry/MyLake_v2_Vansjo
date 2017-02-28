function [TP_obs,TP_mod,chl_obs,chl_mod, mod_all, input_all, INCA_QC, MyLake_results, Sediment_results] = fn_MyL_application(m_start,m_stop, K_sediments, K_lake, use_INCA, run_INCA, run_ID, clim_ID)
global sed_par_file lake_par_file Eevapor
% This is the main MyLake application configuration file. INCA is a switch
% It is made to run a after the parameter are set by Set_Prior

Eevapor=0;    disp('init ...');

calibration_k_values = [(1:50)',cell2mat(K_sediments(:,1)) ]; % writing sediments parameters file

%% generates unique files

sed_par_file = tempname;
lake_par_file = tempname;

dlmwrite(sed_par_file, calibration_k_values,'delimiter','\t');

%% writing lake parameter file
f = fopen('IO/vansjo_para.txt');
garbage = fgetl(f); % file get line
garbage = fgetl(f); % file get line
data_lake = textscan(f, '%s%f%f%f%s', 63, 'Delimiter', '\t');
fclose(f); % the parameter line (xx,1) + 2 lines gives the location of the paramter in the input txt file.
% array position + 2 = input file line

data_lake{1, 2}(23,1) = K_lake{1}; % I_scDOC
data_lake{1, 2}(46,1) = K_lake{2}; % qy_DOC
data_lake{1, 2}(45,1) = K_lake{3}; % oc_DOC
data_lake{1, 2}(47,1) = K_lake{4}; % k_bod
data_lake{1, 2}(49,1) = K_lake{5}; % theta_bod
data_lake{1, 2}(50,1) = K_lake{6}; % theta_bod_ice
data_lake{1, 2}(17,1) = K_lake{7};
data_lake{1, 2}(31,1) = K_lake{8};
data_lake{1, 2}(32,1) = K_lake{9};
data_lake{1, 2}(53,1) = K_lake{10}; % BOD7  on line 55 of the txt file
data_lake{1, 2}(60,1) = K_lake{11}; % I_scO
data_lake{1, 2}(61,1) = K_lake{12}; % Q10
data_lake{1, 2}(62,1) = K_lake{13}; % wc_factor
data_lake{1, 2}(63,1) = K_lake{14}; % T_ref

fid=fopen(lake_par_file,'wt');
fprintf(fid,'\n\n');
fclose(fid);
dlmwrite(lake_par_file, [[1:63]',data_lake{2},data_lake{3},data_lake{4},(1:63)'],'delimiter','\t','-append'); % 1:63 is the lenght of the parameter file.

%% Specific MyLake application

warning('off', 'all')
lake='Vansjo';
year=1983;
dt = 1.0;

%% If you want MyLake to read INCA outputs use_INCA
% (1) Reading exisitng INCA outputs to prepare MyLake input
% (0) You are running only MyLake, inputs alrady exist in the IO folder
if use_INCA == 1;
   [store_INCAP_input,vanem_INCAP_input, INCA_QC]=fn_INCA_MyL(run_INCA, run_ID, clim_ID, m_start, m_stop);
else
    INCA_QC = 0;
end


% store_INCAP_input
% vanem_INCAP_input


%# ############ This is Vansj� Storefj ##############

parafile=lake_par_file;
initfile='IO/store_init.txt';


if use_INCA == 0
    inputfile='IO/store_INCAP_input_baseline.txt';
    disp('Using existing input')
elseif use_INCA == 1
    inputfile = store_INCAP_input; % setting use_INCA to 2 will look for store_INCAP_input
    disp('Using INCA output')
elseif ischar(use_INCA);
    inputfile=use_INCA;
    disp('Using response surfaces array')
end


%% 24/9/2014 This is MyLake DOCOMO commands ...  disabled for now
%[In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz, In_TPz_sed,In_Chlz_sed,In_O2z,In_FIM,Ice0,Wt,Inflw,...
% Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
%    = modelinputs_v12_1b(m_start,m_stop, initfile, 'duh', inputfile, 'duh', parafile, 'duh', dt);
%
%In_DICz = In_DOCz;%dummy DIC;
%
%  [zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,DICzt,CO2zt,O2zt,O2_sat_relt,O2_sat_abst,BODzt,Qzt_sed,lambdazt,...
%        P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,surfaceflux,oxygenflux,CO2_eqt,~,O2_eqt,K0_O2t,CO2_ppmt,dO2Chlt,dO2BODt,dO2SODt,testi1t,testi2t,testi3t, sediments_data_basin2] = ...
%      solvemodel_v12_1b_ut(m_start,m_stop,initfile,'lake', inputfile,'timeseries', parafile,'lake', In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz, In_Chlz,In_DOCz,In_DICz,In_O2z,In_TPz_sed,In_Chlz_sed,In_FIM, ...
%                        Ice0,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names, Bio_par,Bio_par_range,Bio_par_names);
%##################################################################

% note: I removed the DIC/O2 bits here ... take them again from Langtjern
% app when migrating to Mylake DOCOMO
Deposition = 0;

disp('Storefjorden ...')

[In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DICz,In_DOCz,In_TPz_sed,In_Chlz_sed,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_FIM,Ice0,Wt,Inflw,...
    Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names] ...
    = modelinputs_v2(m_start,m_stop, initfile, 'lake', inputfile, 'timeseries', parafile, 'lake', dt);


[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,DICzt,CO2zt,O2zt,NO3zt,NH4zt,SO4zt,HSzt,H2Szt,Fe2zt,Ca2zt,pHzt,CH4zt,Fe3zt,Al3zt,SiO4zt,SiO2zt,diatomzt,O2_sat_relt,O2_sat_abst,BODzt,Qzt_sed,lambdazt,...
    P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,surfaceflux,O2fluxt,CO2_eqt,K0t,O2_eqt,K0_O2t,...
    CO2_ppmt,dO2Chlt,dO2BODt,testi1t,testi2t,testi3t,...
    MyLake_results_basin1, sediment_data_basin1] ...
    = solvemodel_v2(m_start,m_stop,initfile,'lake',inputfile,'timeseries', parafile,'lake',In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_TPz_sed,In_Chlz_sed,In_FIM, ...
    Ice0,Wt,Inflw,Phys_par,Phys_par_range,Phys_par_names, Bio_par,Bio_par_range,Bio_par_names, Deposition);

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


% %# ############ This is Vansj� Vanemfj. ##############
% if isnumeric(use_INCA) % to avoid running two basins in case of RS analysis. 
    
%     if use_INCA == 0
%         land_to_vanem = 'IO/vanem_INCAP_input_baseline.txt'
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

%% returns observed and simulated

%% Load and match obs and sim ... commented for now. Later they need to match sim

tlims=[datenum(m_start):datenum(m_stop)]'; % creating a vector of datums
zinx=find(zz<4); %depth layer considered first 4m because of sampling device

if ischar(use_INCA) % output for response surfaces run
    
    temp=(Pzt(zinx,:)+PPzt(zinx,:)+DOPzt(zinx,:)+Chlzt(zinx,:))'; % total P is computed
    TP_mod = [tlims mean(temp,2)]; % this is the mean of depth-averaged time series
    
    temp=(Chlzt(zinx,:)+Czt(zinx,:))'; % this is the mean of depth-averaged time series
    chl_mod=[tlims mean(temp,2)];
    
    
    TP_obs = [];
    chl_obs = [];
    
end

% cd ..;cd ..;cd ..;cd ..;

if isnumeric(use_INCA)
    

    cd 'obs'; cd 'vanem_obs';
    load 'TOTP.dat' % these are just C&P of vanem ...
    temp=(Pzt(zinx,:)+PPzt(zinx,:)+DOPzt(zinx,:)+Chlzt(zinx,:))'; % total P is computed
    TP_mod_all = [tlims mean(temp,2)]; % this is the mean of depth-averaged time series
    TP_obs = TOTP; % retreiving the observations
    [v,loc_obs,loc_sim] = intersect(TP_obs(:,1), TP_mod_all(:,1)); % returns the datum, and the index for both obs and sim
    MatchedData = [v TP_obs(loc_obs,2) TP_mod_all(loc_sim,2)];% and create subset of data with elements= Time, Observed, Simulated
    TP_obs = MatchedData (:,2); % fn output
    TP_mod = MatchedData (:,3); % fn output
    clear MatchedData
    
    % file_path = ['obs' filesep 'vanem_obs' filesep 'Cha.dat']
    load 'Cha.dat' % these are just C&P of vanem ...
    temp=(Chlzt(zinx,:)+Czt(zinx,:))'; % this is the mean of depth-averaged time series
    chl_mod_all=[tlims mean(temp,2)];
    chl_obs = Cha;
    [v,loc_obs,loc_sim] = intersect(chl_obs(:,1), chl_mod_all(:,1)); % returns the datum, and the index for both obs and sim
    MatchedData = [v chl_obs(loc_obs,2) chl_mod_all(loc_sim,2)];% and create subset of data with elements= Time, Observed, Simulated
    chl_obs = MatchedData (:,2); % fn output
    chl_mod = MatchedData (:,3); % fn output
    clear MatchedData
    
    % file_path = ['obs' filesep 'vanem_obs' filesep 'PO4.dat']
    load 'PO4.dat' % these are just C&P of vanem ...
    temp=Pzt(zinx,:)'; % this is the mean of depth-averaged time series
    PO4_mod_all=[tlims mean(temp,2)];
    PO4_obs = PO4;
    [v,loc_obs,loc_sim] = intersect(PO4_obs(:,1), PO4_mod_all(:,1)); % returns the datum, and the index for both obs and sim
    MatchedData = [v PO4_obs(loc_obs,2) PO4_mod_all(loc_sim,2)];% and create subset of data with elements= Time, Observed, Simulated
    PO4_obs = MatchedData (:,2); % fn output
    PO4_mod = MatchedData (:,3); % fn output
    clear MatchedData
    
    % file_path = ['obs' filesep 'vanem_obs' filesep 'Part.dat']
    load 'Part.dat' % these are just C&P of vanem ...
    temp=PPzt(zinx,:)'; % this is the mean of depth-averaged time series
    Part_mod_all=[tlims mean(temp,2)];
    Part_obs = Part;
    [v,loc_obs,loc_sim] = intersect(Part_obs(:,1), Part_mod_all(:,1)); % returns the datum, and the index for both obs and sim
    MatchedData = [v Part_obs(loc_obs,2) Part_mod_all(loc_sim,2)];% and create subset of data with elements= Time, Observed, Simulated
    Part_obs = MatchedData (:,2); % fn output
    Part_mod = MatchedData (:,3); % fn output
    clear MatchedData
    
    temp=Tzt(zinx,:)'; % this is the mean of depth-averaged time series
    Temp_mod_all=[tlims mean(temp,2)];
    
end    

% MyLake_results = [MyLake_results_basin1, MyLake_results_basin2];
MyLake_results = [MyLake_results_basin1];
% Sediment_results = [sediment_data_basin1, sediment_data_basin2];
Sediment_results = [sediment_data_basin1];

mod_all =  0; %[TP_mod_all, chl_mod_all(:,2), PO4_mod_all(:,2), Part_mod_all(:,2) Temp_mod_all(:,2) Q_lake]; % P speciation epilimnion 2m
input_all = 0; %[Wt,Inflw]; % weather and inflows 

%% cleaning
fclose('all');
% delete (sed_par_file)
% delete (lake_par_file)

end