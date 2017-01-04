clearvars x y
low_cutoff = 1;
x = big_results{1,2}(:,3);
x(x<low_cutoff)= NaN;
col = [2,5,6,7,8,9,10,11,12,13,14]; % those colums with 27759 lines
for i = 2:11
    y = big_results{1,col(i)}(:,3);
    y(y<low_cutoff) = NaN;
    x = [x, y];
end
boxplot(x)