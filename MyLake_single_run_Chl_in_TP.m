% for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[2000, 1, 1]; %
m_stop=[2013, 12, 31]; %

save_initial_conditions = false; % save final concentrations as initial for the next run


[lake_params, sediment_params] = load_params();
name_of_scenario = 'IO/Scenarios/T_only_full_scen_base_historical_20y.txt'

% Niva results RMSD = 130 =======================================
file_name = 'IO/chl_in_tp.mat'
lake_params{47} = 58.3842e-003; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = 128.2949e-003; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} = 1.4988e+000; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = 1.6945e+000; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} = 208.3324e-003; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = 201.6135e-003; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = 1.2687e+000; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = 1.6142e+000; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} = 31.3665e-003; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} = 14.4699e-006; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = 30.5827e-006; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = 37.9560e-003; % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = 34.7141e-003; % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
sediment_params{52} = 21.5114e+000; % 65.1237e+000   %    accel
lake_params{24} = 373.1228e-003; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)


% Trials:
lake_params{24} = 1; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
lake_params{47} = 0.07; % 47     settling velocity for Chl1 a (m day-1)
lake_params{46} = 0.05; % 46  settling velocity for S (m day-1)
lake_params{56} = 0.07; % 56    Settling velocity for Chl2 a (m day-1)
sediment_params{52} = 100; % 65.1237e+000   %    accel
% =====================================================================================================================


% Latest Niva calibration with RMSD 133: % ====================================================
% file_name = 'IO/calib_17_11_01_Niva_calib_rmsd_3x_1x_133.mat'
% x = [0.0691969351707719; 0.106384525900796; 1.33423213718069; 1.81948427216433; 0.0992004077859556; 0.211874114480601; 1.08692651868250; 0.662212213886266; 0.0426104938223779; 2.28258392431412e-05; 7.44071863750726e-05; 0.0432478460870351; 0.0232570145845245; 22.0651718079512; 0.735236639414830];

% Niva results "-RMSD*(R^2 - 1) = 50" % ====================================================
% file_name = 'IO/niva_RMSDxR_2017_11_10.mat'
% x = [0.274263823391964; 0.197850527680004; 1.11383214890828; 1.48285004902205; 0.0646703692387376; 0.113598553694002; 1.38507986430430; 0.846610778723349; 0.0522258868020319; 5.69697106295151e-05; 5.10451900159367e-05; 0.0197061167375151; 0.0424292349113952; 17.5765287697599; 0.680178371230502];

% Ecomac-2 results "-RMSD*(R^2 - 1) = 42" % ====================================================
% file_name = 'IO/test.mat'
% file_name = 'IO/ecomac_2_RMSDxR_2017_11_10.mat'
% x = [0.0509799553636229; 0.112442493655332; 1.28362330332034; 1.36809570914136; 0.0501864460452149; 0.108370399688849; 1.46968213451174; 1.53127867204317; 0.0473597449781419; 2.37096381881603e-05; 3.26782527521984e-05; 0.0449499936004535; 0.0403445049032881; 20.3036606042011; 0.622744575375964];


% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
% =====================================================================================================================


% Ecomac results "-(R^2 - 1) = 4.8027" % ====================================================
% file_name = 'IO/test.mat'
% x = [0.132606917590099, 0.225737487729358, 1.19511818051128, 1.38191729348643, 0.0814948456463502, 0.168250382189148, 1.33924032035532, 0.352470118359444, 0.0638216165399783, 3.92379970293085e-05, 3.29745381818762e-05, 0.0365888738061591, 0.0447918251355218, 20.6787702244692, 0.875859748737663];

% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
% % =====================================================================================================================

% Niva sediment cores = r^2*RMSD, res=962.4  % ====================================================
% file_name = 'IO/test.mat'
% x = [0.174148347056693; 0.237154770484129; 1.25253945867738; 1.97462907422698; 0.0500000000000000; 0.239621978486234; 1.43696261225620; 0.961686642123568; 0.0949427992591365; 6.58225239993542e-05; 2.95856455661290e-05; 0.0102276064570374; 0.00797357868464694; 3.69374890450473; 0.916622502253878; 0.322613371643477; 0.0999920378328245; 0.0245579518367847; 0.0893829812780243; 0.0549224128647664; 0.00132750427627326; 59.7003093123739];

% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
% sediment_params{1} = x(16);  %   'k_Chl',                 %        % 1
% sediment_params{2} = x(17);  %  'k_POP',                 %        % 1
% sediment_params{3} = x(18);  % 'k_POC',                  %        % 0.01
% sediment_params{4} = x(19);  %  'k_DOP',                 %        % 1
% sediment_params{5} = x(20);  % 'k_DOC',                  %        % 1
% sediment_params{23} = x(21);  %     'k_pdesorb_a',         %
% sediment_params{24} = x(22);  %     'k_pdesorb_b',         %
% % =====================================================================================================================


% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=960.4  % ====================================================
% file_name = 'IO/test.mat'
% x = [0.292867131719439; 0.152596930300795; 1.35441691707732; 0.404760282391585; 0.490917403675707; 0.258956557269084; 1.33413227748294; 0.622799283023575; 0.416602967265499; 1.07101519222963e-05; 2.32178478578498e-05; 0.0339306144766581; 0.0299736636364725; 28.2162698311715; 0.935456804191259; 0.395784453077929; 0.0121700676481273; 0.0783899484741691; 0.0834067833948156; 0.0269273099451162; 0.547585590824826; 30.4046816190090; 65.1459885997294; 84.4123799839751; 1.63199703067483; 9.65075155669598; 94.5194422916262; 59.8123208815519; 11.7793522493300; 56.8609132402850; 21.4807095781679; 37.1702379722574; 1.03583807872418];

% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)

% % new added for cores
% sediment_params{1} = x(16);  %   'k_Chl',                 %        % 1
% sediment_params{2} = x(17);  %  'k_POP',                 %        % 1
% sediment_params{3} = x(18);  % 'k_POC',                  %        % 0.01
% sediment_params{4} = x(19);  %  'k_DOP',                 %        % 1
% sediment_params{5} = x(20);  % 'k_DOC',                  %        % 1
% sediment_params{23} = x(21);  %     'k_pdesorb_a',         %
% sediment_params{24} = x(22);  %     'k_pdesorb_b',         %

% % for cores too (scaling unknown inputs):
% lake_params{18} = x(23);%    scaling factor for inflow concentration of C (-)
% lake_params{19} = x(24);%    scaling factor for inflow concentration of POC (-)
% lake_params{20} = x(25);%    scaling factor for inflow concentration of total P (-)
% lake_params{21} = x(26);%    scaling factor for inflow concentration of diss. organic P (-)
% lake_params{22} = x(27);%    scaling factor for inflow concentration of Chl a (-)
% lake_params{23} = x(28);%    scaling factor for inflow concentration of DOC  (-)
% lake_params{25} = x(29);%    Scaling factor for inflow concentration of O2 (-)
% lake_params{27} = x(30);%    Scaling factor for inflow concentration of NO3 (-)
% lake_params{34} = x(31);%    Scaling factor for inflow concentration of Fe3 (-)
% lake_params{35} = x(32);%    Scaling factor for inflow concentration of Al3 (-)
% lake_params{37} = x(33);%    Scaling factor for inflow concentration of CaCO3 (-)
% % =====================================================================================================================


% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=415  % ====================================================
% file_name = 'IO/ecomac-2_.mat'
% x = [0.0954507159077614; 0.200583020651291; 1.44777416925993; 0.444680762956579; 0.0500000000000000; 0.209087437993630; 1.26668516980844; 1.20000000000000; 0.0765004968774690; 1.00000000000000e-05; 1.00000000000000e-05; 0.0108503379659131; 0.0208594885646256; 2; 1; 0.111452954004814; 0.0756702079558015; 0.0633637582049403; 0.100000000000000; 0.0957726460081909; 1.96425357122318; 19.1927004493911; 63.3645672158996; 82.0717029692034; 1.58775222593736; 1; 71.8011131097223; 29.5660628177359; 1.06765260957026; 0.206574107548012; 1; 3.94332721070197; 1.69575727538507];

% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)

% % new added for cores
% sediment_params{1} = x(16);  %   'k_Chl',                 %        % 1
% sediment_params{2} = x(17);  %  'k_POP',                 %        % 1
% sediment_params{3} = x(18);  % 'k_POC',                  %        % 0.01
% sediment_params{4} = x(19);  %  'k_DOP',                 %        % 1
% sediment_params{5} = x(20);  % 'k_DOC',                  %        % 1
% sediment_params{23} = x(21);  %     'k_pdesorb_a',         %
% sediment_params{24} = x(22);  %     'k_pdesorb_b',         %

% % for cores too (scaling unknown inputs):
% lake_params{18} = x(23);%    scaling factor for inflow concentration of C (-)
% lake_params{19} = x(24);%    scaling factor for inflow concentration of POC (-)
% lake_params{20} = x(25);%    scaling factor for inflow concentration of total P (-)
% lake_params{21} = x(26);%    scaling factor for inflow concentration of diss. organic P (-)
% lake_params{22} = x(27);%    scaling factor for inflow concentration of Chl a (-)
% lake_params{23} = x(28);%    scaling factor for inflow concentration of DOC  (-)
% lake_params{25} = x(29);%    Scaling factor for inflow concentration of O2 (-)
% lake_params{27} = x(30);%    Scaling factor for inflow concentration of NO3 (-)
% lake_params{34} = x(31);%    Scaling factor for inflow concentration of Fe3 (-)
% lake_params{35} = x(32);%    Scaling factor for inflow concentration of Al3 (-)
% lake_params{37} = x(33);%    Scaling factor for inflow concentration of CaCO3 (-)
% % % =====================================================================================================================

% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=213  % ====================================================
% file_name = 'IO/niva_cores_inputs_2017_12_11.mat'
% x = [0.0500000000000000; 0.100000000000000; 1; 1.20000000000000; 0.0500000000000000; 0.211177491633875; 1; 2; 1; 0.000100000000000000; 1.00000000000000e-05; 0.00500000000000000; 0.0394039349761114; 35.7967175968888; 0; 0.0100000000000000; 0.0148658880420987; 0.00100000000000000; 0.0564425412638617; 0.100000000000000; 0.152903147742124; 1.00100000000000; 74.8258807476834; 35.9440655491842; 0; 2; 48.7137529836835; 0.962954178292867; 1.09937378894996; 0.861282785490013; 11.3758119915549; 0; 0.250000000000000];

% lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = x(14); % 65.1237e+000   %    accel
% lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)

% % new added for cores
% sediment_params{1} = x(16);  %   'k_Chl',                 %        % 1
% sediment_params{2} = x(17);  %  'k_POP',                 %        % 1
% sediment_params{3} = x(18);  % 'k_POC',                  %        % 0.01
% sediment_params{4} = x(19);  %  'k_DOP',                 %        % 1
% sediment_params{5} = x(20);  % 'k_DOC',                  %        % 1
% sediment_params{23} = x(21);  %     'k_pdesorb_a',         %
% sediment_params{24} = x(22);  %     'k_pdesorb_b',         %

% % for cores too (scaling unknown inputs):
% lake_params{18} = x(23);%    scaling factor for inflow concentration of C (-)
% lake_params{19} = x(24);%    scaling factor for inflow concentration of POC (-)
% lake_params{20} = x(25);%    scaling factor for inflow concentration of total P (-)
% lake_params{21} = x(26);%    scaling factor for inflow concentration of diss. organic P (-)
% lake_params{22} = x(27);%    scaling factor for inflow concentration of Chl a (-)
% lake_params{23} = x(28);%    scaling factor for inflow concentration of DOC  (-)
% lake_params{25} = x(29);%    Scaling factor for inflow concentration of O2 (-)
% lake_params{27} = x(30);%    Scaling factor for inflow concentration of NO3 (-)
% lake_params{34} = x(31);%    Scaling factor for inflow concentration of Fe3 (-)
% lake_params{35} = x(32);%    Scaling factor for inflow concentration of Al3 (-)
% lake_params{37} = x(33);%    Scaling factor for inflow concentration of CaCO3 (-)
% % % % =====================================================================================================================

% Chl in TP results RMSD = 136 =======================================
% file_name = 'IO/chl_in_tp.mat'
lake_params{47} = 4.1304e-01;
lake_params{49} = 1.0431e-01;
lake_params{50} = 1.4133e+00;
lake_params{53} = 1.0170e+00;
lake_params{56} = 7.2653e-02;
lake_params{57} = 1.0026e-01;
lake_params{58} = 1.3597e+00;
lake_params{59} = 1.2248e+00;
lake_params{46} = 7.6604e-01;
lake_params{10} = 5.1714e-05;
lake_params{54} = 1.6050e-05;
lake_params{12} = 5.0137e-03;
lake_params{55} = 4.4973e-02;
sediment_params{52} = 2.8081e+01;
lake_params{24} = 7.8327e-01;


% try
run_ID = 0;
clim_ID = 0;
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, name_of_scenario, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim % runs the model and outputs obs and sim


disp('Saving results...')
save(file_name, 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));

if is_metrics == true

    load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

    depths = [5;10;15;20;25;30;35;40];
    rmsd_O2 = 0;


    for i=1:size(depths,1)
        d = depths(i);
        zinx=find(MyLake_results.basin1.z == d);
        O2_measured = res.T(res.depth1 == d);
        day_measured = res.date(res.depth1 == d);
        day_measured = day_measured(~isnan(O2_measured));
        O2_measured = O2_measured(~isnan(O2_measured));

        O2_mod = MyLake_results.basin1.concentrations.O2(zinx,:)'/1000;
        [T_date,loc_sim, loc_obs] = intersect(MyLake_results.basin1.days, day_measured);

        % rmsd_O2 = rmsd_O2 + RMSE(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        rmsd_O2 = rmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    end

    zinx=find(MyLake_results.basin1.z<4);
    TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.DOP(zinx,:) + MyLake_results.basin1.concentrations.POP(zinx,:))', 2);
    Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
    P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
    POP_mod = mean((MyLake_results.basin1.concentrations.POP(zinx,:) + MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

    load 'obs/store_obs/TOTP.dat' % measured
    % load 'obs/store_obs/Cha.dat' % measured
    load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
    load 'obs/store_obs/PO4.dat' % measured
    load 'obs/store_obs/Part.dat' % measured


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
    rmsd_TOTP = rmsd(TP_mod(loc_sim, 1), TOTP(loc_obs, 2))


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
    rmsd_Chl = rmsd(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    rmsd_PO4 = rmsd(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    rmsd_PP = rmsd(POP_mod(loc_sim, 1), Part(loc_obs, 2));

    disp('RMSD 3xRMSE(P)+RMSE(O2):')
    disp(sum([3*rmsd_TOTP, 3*rmsd_Chl, 3*rmsd_PO4, 3*rmsd_PP, rmsd_O2]))
    disp('RMSD = RMSE(P)+RMSE(O2):')
    disp(sum([rmsd_TOTP, rmsd_Chl, rmsd_PO4, rmsd_PP, rmsd_O2]))
end


toc



% end
%
