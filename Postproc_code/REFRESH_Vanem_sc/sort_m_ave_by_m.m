%% This will make monthly-averaged box-plots

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
   
    for m = 1:3 % mamagenemt scenario loop
        for c = 1:4 % climate scenario loop
            % now we are looking at a 8412 x 1 vector, we want monthly averages
            % by year, that is a 8412/365 x 12 matrix. We could reshape  like so:
            % X = reshape(x,365,Ndays/365), where Ndays = length(temp{m,c})
            % but we dont have an integer number of years ... so we will look
            % ...
            
            % now we will contacenate date_index and our lovely column row vectors
            x_average_vs_t{m,c}=horzcat(date_index,temp{m,c});
            
            %% Loop throug months
            
            %finding starting and ending dates
            starty = x_average_vs_t{m,c}(1,1);
            endy = x_average_vs_t{m,c}(length(temp{1,1}),1);
            startm = x_average_vs_t{m,c}(1,2);
            endm = x_average_vs_t{m,c}(length(temp{1,1}),2);
            startd = x_average_vs_t{m,c}(1,3);
            endd = x_average_vs_t{m,c}(length(temp{1,1}),3);
            
            monthly_average = zeros(1,1); % our average will be stored here
            
            for year = 1:(endy-starty)+1
                current_year = starty+(year-1);
                fprintf('scenario %d-%d, year loop %d\n', m, c, current_year);
                for month = 1:12
                    cumulative_concentration = 0;
                    m_lenght = 0;
                    for datum = 1:8214; % wrong: should be universal
                        if x_average_vs_t{m,c}(datum,1) == current_year
                            if x_average_vs_t{m,c}(datum,2)== month
                                cumulative_concentration = cumulative_concentration + x_average_vs_t{m,c}(datum,4);
                                m_lenght = m_lenght + 1;
                            end
                        end
                    end
                    monthly_average(year,month) = cumulative_concentration/m_lenght;
                end
                
            end
            
            agg_month_ave{m,c}= monthly_average;
            
        end
    end
    
    fprintf('saving arrays %d of 3 ...\n',variable)
    
       
    if variable == 1
        tp_x_average_vs_t = x_average_vs_t;
        tp_agg_month_ave = agg_month_ave;
    elseif variable == 2
        DIP_agg_month_ave = agg_month_ave;
        DIP_x_average_vs_t = x_average_vs_t;
    else
        POP_agg_month_ave = agg_month_ave;
        POP_x_average_vs_t = x_average_vs_t;
    end
    
end
clear x_average_vs_t agg_month_ave monthly_average date_index temp c cumulative_concentration datum endd endm endy m m_lenght startd startm starty variable year current_year month
clear DIPdepths POPdepths totalPdepths volumeflowout 
save agg_month_ave.mat
fprintf('done.\n')
