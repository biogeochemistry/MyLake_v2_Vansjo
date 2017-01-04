function agg = fn_agg_t_depth_averaged(var_name, d_up, d_bottom)
% This function takes the result array for a specified variable and aggerates
% into a single date (year, month, day) vs concentration matrix. For not it is specific to TP

cd('C:\skydrive\PROJECT\O-10020 REFRESH\SV_output\Vanem_sc');
global totalPdepths DIPdepths POPdepths;

conc = zeros(8214,19); %concentration matrix
pset_temp = cell(45,1); % vector array of parameter sets
working_array = cell(3,4); % matrix of scenarios
%% load weighting vector. modify to change param. weight

load('param_w.mat', 'param_w')

%% averages conconcentration over depth, and weighted parameter set


for m = 1:3 % mamagenemt scenario loop
    for c = 1:4 % climate scenario loop
        
%% averaging for desired depth for each individual parameter sets        
        for p = 1:45 % parameter set loop
            % concentration vs depth matrix becomes concentration col_vec
            if var_name == 1
                d_mat=totalPdepths{m,c}{p,1}(:,[d_up:d_bottom]); %gathers desired water depths
            elseif var_name == 2
                d_mat=DIPdepths{m,c}{p,1}(:,[d_up:d_bottom]);
            else
                d_mat=POPdepths{m,c}{p,1}(:,[d_up:d_bottom]);
            end
            
            temp1=mean(d_mat'); % concentration col_vec
            pset_temp{p}=temp1'; % concentration col_vev for each pset 
            
            
        end
        
        temp2 = {1:length(param_w)};
%% weigned average of parameter set, 10 and 90 pecentiles
        
% multiply pset values by weight
        for i = 1: length(param_w) %index param vector
            temp2{i}=(pset_temp{i}*param_w(i)); % temp2{} containes conc vs time weighted for psets
        end
        
        %pset_temp{} is an array of col_vec containing conc vs time
        %we make a matrix with pset_temp
        pset_mat = (horzcat(pset_temp{:,:}));
        
        %now we can do stat on pset_mat
        %the quartile function we made only works with column ...
        for datum = 1:8214  %generalize
        q_stat = fn_quartile((pset_mat(datum,:)'), 15,85);
        q_stat_col(datum,1) = q_stat(1);
        q_stat_col(datum,2) = q_stat(2);
        end
        
        weighted_ave = (sum(cell2mat(temp2)')')/60; % creates col_vec with weighted average
        working_array{m,c}(:,1)= weighted_ave; % weighted conc col_vc placed in {c,m} array
        working_array{m,c}(:,2)= q_stat_col(:,1);
        working_array{m,c}(:,3)= q_stat_col(:,2);
        
   %     working_array{m,c}(:,2)= 
        
        
    end
end

agg = working_array;