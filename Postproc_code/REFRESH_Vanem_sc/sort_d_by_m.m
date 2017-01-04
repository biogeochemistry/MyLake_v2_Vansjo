%% This will prepare box-plots for each months, with daily information in them

%% Asks for the depth interval over which to calculate depth-averaged concentrations
% we want a time vs concentrations result ... not depth profile
concacenated_with_dates = cell(3,4); % we need this one to concacenate date and data
%fprintf('loading vardepth.mat ... \n')
%load vardepth.mat
%fprintf('done \n')
load('date_index.mat')
global totalPdepths DIPdepths POPdepths;

%% We will do a loop over the three parameters
for variable = 1:3
    % Calls the aggegator, returns {c,m} array whose elements are depth averaged, p-weighter conc. time-series
    temp = fn_agg_t_depth_averaged(variable, 1, 4); % variable name, top boundary, botton boundary for average
    % temp will contain an {m,c} array with concentrations at each timeteps
    for m = 1:3 % mamagenemt scenario loop
        for c = 1:4 % climate scenario loop
            % now we are looking at a 8412 x 1 vector, we want monthly averages
            % by year, that is a 8412/365 x 12 matrix. We could reshape  like so:
            % X = reshape(x,365,Ndays/365), where Ndays = length(temp{m,c})
            % but we dont have an integer number of years ... so we will look
            % ...
            
            % now we will contacenate date_index and our col_vec made by
            % fn_agg_t_depth_average
            
            x_average_vs_t{m,c}=horzcat(date_index,temp{m,c});
            
            %% Loop throug months
            
            %finding starting and ending dates
            starty = x_average_vs_t{m,c}(1,1);
            endy = x_average_vs_t{m,c}(length(temp{1,1}),1);
            startm = x_average_vs_t{m,c}(1,2);
            endm = x_average_vs_t{m,c}(length(temp{1,1}),2);
            startd = x_average_vs_t{m,c}(1,3);
            endd = x_average_vs_t{m,c}(length(temp{1,1}),3);
            
            for month = 1:12
                day_token = 1;
                for datum = 1:8214;
                    if x_average_vs_t{m,c}(datum,2)== month
                        daily_value(month,day_token) = x_average_vs_t{m,c}(datum,4);
                        day_token = day_token + 1;
                    end
                end
            end
            agg_ave{m,c}= daily_value;
        end
    end
    fprintf('saving arrays %d of 3 ...\n',variable)
    % saves the arrays to the right simulated variable
    if variable == 1
        tp_x_average_vs_t = x_average_vs_t;
        tp_day_by_month = agg_ave;
    elseif variable == 2
        DIP_day_by_month = agg_ave;
        DIP_x_average_vs_t = x_average_vs_t;
    else
        POP_day_by_month = agg_ave;
        POP_x_average_vs_t = x_average_vs_t;
    end
end
clear x_average_vs_t agg_ave daily_value date_index temp c cumulative_concentration datum endd endm endy m m_lenght startd startm starty variable year current_year month
%clear DIPdepths POPdepths totalPdepths volumeflowout
%save d_by_m.mat

 fprintf('inverting, removing zeros ...')
%% prepare for cut&paste to sigmaplot 
% one box-plot for each month, contains day values

token = 1;
for c = 1:4 % mamagenemt scenario loop
    for m = 1:3 % climate scenario loop
        
        temp = tp_day_by_month{m,c}';
        sigma_mat_m_tp(:,token:token+11) = temp(:, :);     
        
        temp = POP_day_by_month{m,c}';
        sigma_mat_m_POP(:,token:token+11) = temp(:, :);     
        
        temp = DIP_day_by_month{m,c}';
        sigma_mat_m_DIP(:,token:token+11) = temp(:, :);     
                
        token = token + 12;
    end
end

sigma_mat_m_PIP = sigma_mat_m_tp - sigma_mat_m_DIP;

%removing zeroes from box-plots 

sigma_mat_m_tp(sigma_mat_m_tp == 0) = NaN; %note the elegent way of replacing any given values in a big matrix 
sigma_mat_m_DIP(sigma_mat_m_DIP == 0) = NaN;  
sigma_mat_m_POP(sigma_mat_m_POP == 0) = NaN; 
sigma_mat_m_PIP(sigma_mat_m_PIP == 0) = NaN; 

fprintf('done.\n')
