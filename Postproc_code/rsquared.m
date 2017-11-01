function res = rsquared(s, o)
    o_mean = mean(o);
    se = squared_error(s, o);
    se_mean = squared_error(o, o_mean);
    res = 1 - (se / se_mean);
end

