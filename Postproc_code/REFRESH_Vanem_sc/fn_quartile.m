function [stat] = quartile(col_vec, low, high)

% define data set
x = col_vec;
Nx = size(x,1);

%% Percentile Calculation Example

% define percent
kpercent = low;
% STEP 1 - rank the data
y = sort(x);
% STEP 2 - find k% (k /100) of the sample size, n.
k = kpercent/100;
result = k*Nx;
% STEP 3 - if this is an integer, add 0.5. If it isn't an integer round up.
[N,D] = rat(k*Nx);
if isequal(D,1),               % k*Nx is an integer, add 0.5
    result = result+0.5;
else                           % round up
    result = round(result);
end

% STEP 4 - Find the number in this position. If your depth ends 
% in 0.5, then take the midpoint between the two numbers.
[T,R] = strtok(num2str(result),'0.5');
if strcmp(R,'.5'),
    Qk = mean(y(result-0.5:result+0.5));
else
    Qk = y(result);
end
q_low = Qk;


% define percent
kpercent = high;
% STEP 1 - rank the data
y = sort(x);
% STEP 2 - find k% (k /100) of the sample size, n.
k = kpercent/100;
result = k*Nx;
% STEP 3 - if this is an integer, add 0.5. If it isn't an integer round up.
[N,D] = rat(k*Nx);
if isequal(D,1),               % k*Nx is an integer, add 0.5
    result = result+0.5;
else                           % round up
    result = round(result);
end

% STEP 4 - Find the number in this position. If your depth ends 
% in 0.5, then take the midpoint between the two numbers.
[T,R] = strtok(num2str(result),'0.5');
if strcmp(R,'.5'),
    Qk = mean(y(result-0.5:result+0.5));
else
    Qk = y(result);
end
q_high = Qk;

stat = [q_low q_high];