% This script reads compressed Mylake Output from amazon clusterand creates a nice MATLAB array with it.

clear
cd('C:\landrive\REFRESH\SV_output\Vanem_sc') %we are in the right folder now

letter = input('Unzip files (y/n)?', 's');

%% Deciding on the variable to work with (one x,t array, one t vector)
% 
varname1 = 'totalPdepths';
varname2 = 'volumeflowout';

varname3 = 'DIPdepths';
varname4 = 'POPdepths';


%% Unzipping if necessary
if letter == 'y'
    
    for p = 1:45 % parameter set
        for c = 1:4 % climate scenario
            for m = 1:3 % mamagenemt scenario
                
                if p < 10 % pad the right amount of zeroes
                    dir_name = sprintf('runs\\00%d-c%d_m%d',p,c,m);
                else
                    dir_name = sprintf('runs\\0%d-c%d_m%d',p,c,m);
                end
                
                tp_file = sprintf('%s\\totalPdepths.csv.bz2',dir_name);
                DIP_file = sprintf('%s\\DIPdepths.csv.bz2',dir_name);
                POP_file = sprintf('%s\\POPdepths.csv.bz2',dir_name);
                vol_file = sprintf('%s\\volumeflowout.csv.bz2',dir_name);
                
                unzip = sprintf('bunzip2 -k -q %s',tp_file);
                status = system(unzip);
                
                unzip = sprintf('bunzip2 -k -q %s',DIP_file);
                status = system(unzip);
                
                unzip = sprintf('bunzip2 -k -q %s',POP_file);
                status = system(unzip);
                
                unzip = sprintf('bunzip2 -k -q %s',vol_file);
                status = system(unzip);
                
            end
        end
    end
    
end


%% Pre-allocating the arrays
% for time vs depth concentrations profiles
conc = zeros(8214,19); %concentration matrix

pset_tp = cell(45,1); % vector array of parameter sets
totalPdepths = cell(3,4); % matrix of scenarios

pset_flow = cell(45,1);
volumeflowout = cell(3,4);

pset_DIP = cell(45,1);  
DIPdepths = cell(3,4); 

pset_POP = cell(45,1); 
POPdepths = cell(3,4); 

%% Key
%- lake_x will be a vector {v} 1x9 for each modeled variable. not now.
%- each element of v is a {c,m} 3x4 array for scenarios
%- each element of {c,m} is {pset} a 1x45 vector for param sets
%- each element of {pset} is conc, a (t,d) 8214x19 matrix containing concentrations in
%   the case of a depth profile, or t, a (8214x1) vector in the case of a time
%   variable such as outflow volume.

for m = 1:3 % mamagenemt scenario loop
    for c = 1:4 % climate scenario loop
        for p = 1:45 % parameter set loop
            %% Prepares the filenames and pads with 0 if necessary
            if p < 10
                dir_name = sprintf('runs\\00%d-c%d_m%d',p,c,m);
            else
                dir_name = sprintf('runs\\0%d-c%d_m%d',p,c,m);
            end
            
            %% Places the 00p-cc_mm\files.csv into temporary mat or vec
            tp_file = sprintf('%s\\%s.csv',dir_name, varname1);
            vol_file = sprintf('%s\\%s.csv',dir_name, varname2);

            DIP_file = sprintf('%s\\%s.csv',dir_name, varname3);
            POP_file = sprintf('%s\\%s.csv',dir_name, varname4);
            
            %because some folder are empty  
            if strncmp(dir_name, 'runs\021-c3_m2',13) == 0 && strncmp(dir_name, 'runs\024-c3_m2',13) == 0;
             %see the import functions ...  
                conc = fn_import_depth_profiles(tp_file, 1, 8214);
                pset_tp{p} = conc; % adding the conc matrix to the pset vector
                clear conc
                
                conc = fn_import_depth_profiles(DIP_file, 1, 8214);
                pset_DIP{p} = conc; % adding the conc matrix to the pset vector
                clear conc
                
                conc = fn_import_depth_profiles(POP_file, 1, 8214);
                pset_POP{p} = conc; % adding the conc matrix to the pset vector
                clear conc
                
                vol = fn_import_time_vec(vol_file, 1, 8214);
                pset_flow{p} = vol; % adding the conc matrix to the pset vector
                clear vol  
                
                fprintf('successfully read %s \n',dir_name)
                             
            end
        end
        
        totalPdepths{m,c}=pset_tp;
        DIPdepths{m,c}=pset_DIP;
        POPdepths{m,c}=pset_POP;
        volumeflowout{m,c}=pset_flow;
       
    end
end
fprintf('done ... now saving arrays\n')
clear c letter dir_name m p pset_DIP pset_tp pset_POP pset_flow tp_file DIP_file POP_file vol_file varname1 varname2 varname3 varname4
fprintf('done.\n')