aggregate_yearly_data = cell(3,4);

for m = 1:3 % mamagenemt scenario loop
    for c = 1:4 % climate scenario loop
   
        aggregate_yearly_data{m,c} = aggregate_monthly_averages{m,c}'; 
    end
end