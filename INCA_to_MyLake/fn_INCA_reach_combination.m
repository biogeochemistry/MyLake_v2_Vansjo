function [New_Reach] = fn_INCA_reach_combination(Reach, r_x, r_y, no_vars)
%INCA_REACH_COMBINATION This combines two reaches for Q, C, and Temperature
%   This scripts first combines flow, then loops from 2 to end-1 column,
%   then combines temperature (assuming its always last)

% Q is in m3/s, stuff is in mg/L, temp is in C ... We'll sum flow
New_Reach{1,1}(:,1) = Reach{r_x,1}{1,1}(:,1) + Reach{r_y,1}{1,1}(:,1);
% Stuff in mg/L% ... We'll calculate fluxes of stuff in mg/s for each reach, average to get new reach flux, then divide by new flow to get new concentration
for i = 2:no_vars-1; 
New_Reach{1,i}(:,1) = (((Reach{r_x,1}{1,1}(:,1) .* Reach{r_x,1}{1,i}(:,1) + (Reach{r_y,1}{1,1}(1,1) .* Reach{r_y,1}{1,i}(:,1)))/2 ) ./ New_Reach{1,1}(:,1));
end
% Stuff in celcius ... we'll do a volume base average  
day_sec = 60*60*24;
new_vol = (Reach{r_x,1}{1,1}(:,1) .* day_sec + Reach{r_y,1}{1,1}(:,1) .* day_sec);
New_Reach{1,7}(:,1) = (Reach{r_x,1}{1,1}(:,1) .* day_sec ./ new_vol .* Reach{r_x,1}{1,7}(:,1)) + Reach{r_y,1}{1,1}(:,1) .* day_sec ./ new_vol .* Reach{r_y,1}{1,7}(:,1);

end

